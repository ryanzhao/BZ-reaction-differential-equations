#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

struct Params {
    int nx = 320;
    int ny = 320;
    double dx = 1.0;
    double dt = 0.0022;
    int steps = 260000;
    int display_every = 2;

    double epsilon = 0.04;
    double f = 1.35;
    double q = 0.002;
    double phi = 0.075;

    double Du = 0.55;
    double Dv = 0.32;
    double Dw = 0.10;

    bool save_frames = false;
    int frame_every = 20;
    std::string output_dir = "output";
    int pixel_scale = 2;

    bool pacemaker = true;
    int pacemaker_period = 260;
    int pacemaker_radius = 10;
    int layered_rings = 3;
    int ring_spacing = 12;
    double pacemaker_orbit = 0.10;
    double pacemaker_omega = 0.012;
    int wavebreak_period = 900;
    double inward_ring_width = 16.0;

    enum class ViewField { U, W, Mix };
    ViewField view_field = ViewField::Mix;

    enum class WaveMode { Outward, Inward, Dual, Spiral };
    WaveMode wave_mode = WaveMode::Spiral;

    bool auto_contrast = true;
    bool dish_render = true;
    double heterogeneity = 0.08;
    double temporal_mod = 0.03;
    int secondary_centers = 10;
    int secondary_start = 2800;
    double secondary_strength = 0.48;
};

struct Field {
    int nx;
    int ny;
    std::vector<double> data;

    Field(int nx_, int ny_, double init = 0.0) : nx(nx_), ny(ny_), data(nx_ * ny_, init) {}

    double& operator()(int i, int j) { return data[i + nx * j]; }
    const double& operator()(int i, int j) const { return data[i + nx * j]; }
};

static inline int clamp_idx(int x, int lo, int hi) {
    return std::max(lo, std::min(hi, x));
}

static double laplacian(const Field& a, int i, int j, double dx) {
    int il = clamp_idx(i - 1, 0, a.nx - 1);
    int ir = clamp_idx(i + 1, 0, a.nx - 1);
    int jb = clamp_idx(j - 1, 0, a.ny - 1);
    int jt = clamp_idx(j + 1, 0, a.ny - 1);
    return (a(il, j) + a(ir, j) + a(i, jb) + a(i, jt) - 4.0 * a(i, j)) / (dx * dx);
}

static double view_value(const Params& p, double u, double w) {
    if (p.view_field == Params::ViewField::U) {
        return u;
    }
    if (p.view_field == Params::ViewField::W) {
        return w;
    }
    return std::clamp(0.7 * u + 0.3 * w, 0.0, 1.0);
}

static bool in_dish(const Params& p, int i, int j, double* r_out = nullptr) {
    const double cx = 0.5 * (p.nx - 1);
    const double cy = 0.5 * (p.ny - 1);
    const double dish_radius = 0.47 * std::min(p.nx, p.ny);
    const double dx = i - cx;
    const double dy = j - cy;
    const double r = std::sqrt(dx * dx + dy * dy);
    if (r_out) {
        *r_out = r / dish_radius;
    }
    return r <= dish_radius;
}

static void map_petri_color(double x, double r_norm, bool in_domain, unsigned char& r, unsigned char& g,
                            unsigned char& b) {
    if (!in_domain) {
        r = 243;
        g = 241;
        b = 239;
        if (r_norm > 0.98 && r_norm < 1.06) {
            r = 228;
            g = 212;
            b = 210;
        }
        return;
    }

    const double wave = std::pow(std::clamp(x, 0.0, 1.0), 0.85);
    const double base_r = 210.0;
    const double base_g = 108.0;
    const double base_b = 122.0;
    const double crest_r = 245.0;
    const double crest_g = 242.0;
    const double crest_b = 250.0;

    const double rr = base_r * (1.0 - wave) + crest_r * wave;
    const double gg = base_g * (1.0 - wave) + crest_g * wave;
    const double bb = base_b * (1.0 - wave) + crest_b * wave;

    r = static_cast<unsigned char>(std::clamp(rr, 0.0, 255.0));
    g = static_cast<unsigned char>(std::clamp(gg, 0.0, 255.0));
    b = static_cast<unsigned char>(std::clamp(bb, 0.0, 255.0));
}

static void build_visual_field(const Params& p, const Field& u, const Field& w, Field& visual) {
    double vmin = 1e9;
    double vmax = -1e9;

    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            if (!in_dish(p, i, j)) {
                visual(i, j) = 0.0;
                continue;
            }
            const double v = view_value(p, u(i, j), w(i, j));
            visual(i, j) = v;
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
    }

    if (!p.auto_contrast || vmax - vmin < 1e-9) {
        return;
    }

    const double inv = 1.0 / (vmax - vmin);
    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            if (in_dish(p, i, j)) {
                visual(i, j) = std::clamp((visual(i, j) - vmin) * inv, 0.0, 1.0);
            }
        }
    }
}

static void write_ppm(const Params& p, const Field& visual, int frame_id, const std::string& out_dir) {
    std::filesystem::create_directories(out_dir);
    std::ostringstream name;
    name << out_dir << "/frame_" << std::setw(6) << std::setfill('0') << frame_id << ".ppm";
    std::ofstream ofs(name.str(), std::ios::binary);
    ofs << "P6\n" << visual.nx << " " << visual.ny << "\n255\n";

    for (int j = 0; j < visual.ny; ++j) {
        for (int i = 0; i < visual.nx; ++i) {
            unsigned char r = 0, g = 0, b = 0;
            double r_norm = 0.0;
            const bool inside = in_dish(p, i, j, &r_norm);
            map_petri_color(visual(i, j), r_norm, inside, r, g, b);
            ofs.write(reinterpret_cast<char*>(&r), 1);
            ofs.write(reinterpret_cast<char*>(&g), 1);
            ofs.write(reinterpret_cast<char*>(&b), 1);
        }
    }
}

static int mask_shift(unsigned long mask) {
    int shift = 0;
    while (mask && ((mask & 1UL) == 0UL)) {
        ++shift;
        mask >>= 1UL;
    }
    return shift;
}

static int mask_bits(unsigned long mask) {
    int bits = 0;
    while (mask) {
        bits += (mask & 1UL) ? 1 : 0;
        mask >>= 1UL;
    }
    return bits;
}

static unsigned long channel_to_masked(unsigned char c, unsigned long mask) {
    if (mask == 0) {
        return 0;
    }
    const int shift = mask_shift(mask);
    const int bits = mask_bits(mask);
    const unsigned long max_val = (1UL << bits) - 1UL;
    const unsigned long scaled = static_cast<unsigned long>((static_cast<double>(c) / 255.0) * max_val + 0.5);
    return (scaled << shift) & mask;
}

class X11Renderer {
   public:
    explicit X11Renderer(const Params& p)
        : params_(p), width_(p.nx * p.pixel_scale), height_(p.ny * p.pixel_scale) {
        display_ = XOpenDisplay(nullptr);
        if (!display_) {
            throw std::runtime_error("无法连接到 X11 显示器。请在图形桌面下运行，或使用 --headless。");
        }

        screen_ = DefaultScreen(display_);
        window_ = XCreateSimpleWindow(display_, RootWindow(display_, screen_), 10, 10, width_, height_, 1,
                                      BlackPixel(display_, screen_), WhitePixel(display_, screen_));
        XStoreName(display_, window_, "BZ reaction (spiral waves)");
        XSelectInput(display_, window_, ExposureMask | KeyPressMask | StructureNotifyMask);
        XMapWindow(display_, window_);

        gc_ = DefaultGC(display_, screen_);
        visual_ = DefaultVisual(display_, screen_);
        depth_ = DefaultDepth(display_, screen_);

        bytes_per_pixel_ = 4;
        image_stride_ = width_ * bytes_per_pixel_;
        image_buffer_.resize(image_stride_ * height_);

        image_ = XCreateImage(display_, visual_, depth_, ZPixmap, 0,
                              reinterpret_cast<char*>(image_buffer_.data()), width_, height_, 32, image_stride_);
        if (!image_) {
            throw std::runtime_error("创建 X11 图像失败。");
        }

        red_mask_ = visual_->red_mask;
        green_mask_ = visual_->green_mask;
        blue_mask_ = visual_->blue_mask;
    }

    ~X11Renderer() {
        if (image_) {
            image_->data = nullptr;
            XDestroyImage(image_);
        }
        if (display_) {
            if (window_) {
                XDestroyWindow(display_, window_);
            }
            XCloseDisplay(display_);
        }
    }

    bool process_events() {
        while (XPending(display_)) {
            XEvent ev;
            XNextEvent(display_, &ev);
            if (ev.type == DestroyNotify || ev.type == KeyPress) {
                return false;
            }
        }
        return true;
    }

    void draw_visual(const Field& visual, int step) {
        for (int j = 0; j < params_.ny; ++j) {
            for (int i = 0; i < params_.nx; ++i) {
                unsigned char r, g, b;
                double r_norm = 0.0;
                const bool inside = in_dish(params_, i, j, &r_norm);
                map_petri_color(visual(i, j), r_norm, inside, r, g, b);
                unsigned long pixel = channel_to_masked(r, red_mask_) | channel_to_masked(g, green_mask_) |
                                      channel_to_masked(b, blue_mask_);

                for (int sy = 0; sy < params_.pixel_scale; ++sy) {
                    const int y = j * params_.pixel_scale + sy;
                    auto* row = reinterpret_cast<std::uint32_t*>(image_buffer_.data() + y * image_stride_);
                    for (int sx = 0; sx < params_.pixel_scale; ++sx) {
                        const int x = i * params_.pixel_scale + sx;
                        row[x] = static_cast<std::uint32_t>(pixel);
                    }
                }
            }
        }

        XPutImage(display_, window_, gc_, image_, 0, 0, 0, 0, width_, height_);
        std::ostringstream title;
        title << "BZ spiral wave step=" << step << " (按任意键退出)";
        XStoreName(display_, window_, title.str().c_str());
        XFlush(display_);
    }

   private:
    Params params_;
    int width_;
    int height_;

    Display* display_ = nullptr;
    int screen_ = 0;
    Window window_ = 0;
    GC gc_ = 0;
    Visual* visual_ = nullptr;
    int depth_ = 24;

    XImage* image_ = nullptr;
    int bytes_per_pixel_ = 4;
    int image_stride_ = 0;
    std::vector<unsigned char> image_buffer_;

    unsigned long red_mask_ = 0;
    unsigned long green_mask_ = 0;
    unsigned long blue_mask_ = 0;
};

static void parse_args(int argc, char** argv, Params& p) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--headless") {
            p.display_every = 0;
        } else if (a == "--save") {
            p.save_frames = true;
        } else if (a.rfind("--steps=", 0) == 0) {
            p.steps = std::stoi(a.substr(8));
        } else if (a.rfind("--dt=", 0) == 0) {
            p.dt = std::stod(a.substr(5));
        } else if (a.rfind("--scale=", 0) == 0) {
            p.pixel_scale = std::max(1, std::stoi(a.substr(8)));
        } else if (a.rfind("--nx=", 0) == 0) {
            p.nx = std::max(32, std::stoi(a.substr(5)));
        } else if (a.rfind("--ny=", 0) == 0) {
            p.ny = std::max(32, std::stoi(a.substr(5)));
        } else if (a.rfind("--frame-every=", 0) == 0) {
            p.frame_every = std::max(1, std::stoi(a.substr(14)));
        } else if (a == "--no-pacemaker") {
            p.pacemaker = false;
        } else if (a.rfind("--pacemaker-period=", 0) == 0) {
            p.pacemaker_period = std::max(20, std::stoi(a.substr(19)));
        } else if (a.rfind("--layers=", 0) == 0) {
            p.layered_rings = std::max(1, std::stoi(a.substr(9)));
        } else if (a.rfind("--ring-spacing=", 0) == 0) {
            p.ring_spacing = std::max(4, std::stoi(a.substr(15)));
        } else if (a.rfind("--orbit=", 0) == 0) {
            p.pacemaker_orbit = std::clamp(std::stod(a.substr(8)), 0.0, 0.45);
        } else if (a.rfind("--omega=", 0) == 0) {
            p.pacemaker_omega = std::clamp(std::stod(a.substr(8)), 0.0, 0.08);
        } else if (a.rfind("--wavebreak-period=", 0) == 0) {
            p.wavebreak_period = std::max(80, std::stoi(a.substr(19)));
        } else if (a.rfind("--inward-width=", 0) == 0) {
            p.inward_ring_width = std::max(2.0, std::stod(a.substr(15)));
        } else if (a.rfind("--heterogeneity=", 0) == 0) {
            p.heterogeneity = std::max(0.0, std::stod(a.substr(16)));
        } else if (a.rfind("--temporal-mod=", 0) == 0) {
            p.temporal_mod = std::max(0.0, std::stod(a.substr(15)));
        } else if (a.rfind("--secondary-centers=", 0) == 0) {
            p.secondary_centers = std::max(0, std::stoi(a.substr(20)));
        } else if (a.rfind("--secondary-start=", 0) == 0) {
            p.secondary_start = std::max(0, std::stoi(a.substr(18)));
        } else if (a.rfind("--secondary-strength=", 0) == 0) {
            p.secondary_strength = std::max(0.0, std::stod(a.substr(21)));
        } else if (a.rfind("--view=", 0) == 0) {
            const std::string v = a.substr(7);
            if (v == "u") {
                p.view_field = Params::ViewField::U;
            } else if (v == "w") {
                p.view_field = Params::ViewField::W;
            } else {
                p.view_field = Params::ViewField::Mix;
            }
        } else if (a.rfind("--wave-mode=", 0) == 0) {
            const std::string v = a.substr(12);
            if (v == "inward") {
                p.wave_mode = Params::WaveMode::Inward;
            } else if (v == "dual") {
                p.wave_mode = Params::WaveMode::Dual;
            } else if (v == "spiral") {
                p.wave_mode = Params::WaveMode::Spiral;
            } else {
                p.wave_mode = Params::WaveMode::Outward;
            }
        } else if (a.rfind("--contrast=", 0) == 0) {
            p.auto_contrast = (a.substr(11) != "fixed");
        }
    }
}

static void apply_source_patch(Field& u, Field& w, int cx, int cy, int radius, double u_peak, double w_peak) {
    for (int j = 0; j < u.ny; ++j) {
        for (int i = 0; i < u.nx; ++i) {
            const int dx = i - cx;
            const int dy = j - cy;
            if (dx * dx + dy * dy <= radius * radius) {
                u(i, j) = std::max(u(i, j), u_peak);
                w(i, j) = std::max(w(i, j), w_peak);
            }
        }
    }
}

static void apply_layered_source(Field& u, Field& w, int cx, int cy, int core_radius, int layers, int spacing,
                                 double u_peak, double w_peak) {
    const double sigma = 2.8;
    for (int j = 0; j < u.ny; ++j) {
        for (int i = 0; i < u.nx; ++i) {
            const double dx = static_cast<double>(i - cx);
            const double dy = static_cast<double>(j - cy);
            const double r = std::sqrt(dx * dx + dy * dy);

            double amp = 0.0;
            if (r <= core_radius) {
                amp = 1.0;
            }
            for (int k = 1; k <= layers; ++k) {
                const double target = core_radius + k * spacing;
                const double ring = std::exp(-((r - target) * (r - target)) / (2.0 * sigma * sigma));
                amp = std::max(amp, ring);
            }

            if (amp > 0.12) {
                u(i, j) = std::max(u(i, j), 0.10 + u_peak * amp);
                w(i, j) = std::max(w(i, j), 0.08 + w_peak * amp);
            }
        }
    }
}

static void apply_boundary_ring_source(const Params& p, Field& u, Field& w, double width, double u_peak,
                                       double w_peak) {
    const double cx = 0.5 * (p.nx - 1);
    const double cy = 0.5 * (p.ny - 1);
    const double dish_radius = 0.47 * std::min(p.nx, p.ny);
    const double inner = std::max(0.0, dish_radius - width);

    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            double rnorm = 0.0;
            if (!in_dish(p, i, j, &rnorm)) {
                continue;
            }
            const double dx = i - cx;
            const double dy = j - cy;
            const double rr = std::sqrt(dx * dx + dy * dy);
            if (rr >= inner && rr <= dish_radius) {
                const double t = (rr - inner) / std::max(1e-9, dish_radius - inner);
                const double amp = std::exp(-6.0 * t * t);
                u(i, j) = std::max(u(i, j), 0.10 + u_peak * amp);
                w(i, j) = std::max(w(i, j), 0.08 + w_peak * amp);
            }
        }
    }
}

static void seed_spiral_break(const Params& p, Field& u, Field& w) {
    const int cx = p.nx / 2;
    const int cy = p.ny / 2;
    const int r0 = p.nx / 8;

    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            const int dx = i - cx;
            const int dy = j - cy;
            const double rr = std::sqrt(static_cast<double>(dx * dx + dy * dy));
            if (rr > r0 - 2 && rr < r0 + 2 && dx > 0) {
                u(i, j) = 1.05;
                w(i, j) = 0.9;
            }
            if (rr > r0 - 1 && rr < r0 + 1 && dy > 0 && dx < 0) {
                // small broken segment -> spiral tips
                u(i, j) = 0.2;
            }
        }
    }
}

static void apply_pacemaker(const Params& p, int step, Field& u, Field& w) {
    if (!p.pacemaker) {
        return;
    }

    const int cx = p.nx / 2;
    const int cy = p.ny / 2;

    if (p.wave_mode == Params::WaveMode::Outward || p.wave_mode == Params::WaveMode::Dual) {
        apply_layered_source(u, w, cx, cy, p.pacemaker_radius, p.layered_rings, p.ring_spacing, 0.85, 0.72);
    }

    if (p.wave_mode == Params::WaveMode::Inward || p.wave_mode == Params::WaveMode::Dual) {
        apply_boundary_ring_source(p, u, w, p.inward_ring_width, 0.93, 0.72);
    }

    if (p.wave_mode == Params::WaveMode::Spiral) {
        const double phase = step * p.pacemaker_omega;
        const double orbit = p.pacemaker_orbit * std::min(p.nx, p.ny);
        const int c1x = static_cast<int>(std::round(cx + orbit * std::cos(phase)));
        const int c1y = static_cast<int>(std::round(cy + orbit * std::sin(phase)));
        apply_layered_source(u, w, c1x, c1y, p.pacemaker_radius, p.layered_rings, p.ring_spacing, 0.82, 0.70);

        if (p.wavebreak_period > 0 && step % p.wavebreak_period == 0) {
            const int bx = static_cast<int>(std::round(cx + 0.12 * p.nx * std::cos(phase + 0.8)));
            const int by = static_cast<int>(std::round(cy + 0.12 * p.ny * std::sin(phase + 0.8)));
            for (int j = std::max(0, by - 6); j <= std::min(p.ny - 1, by + 6); ++j) {
                for (int i = std::max(0, bx - 6); i <= std::min(p.nx - 1, bx + 6); ++i) {
                    const int dx = i - bx;
                    const int dy = j - by;
                    if (dx * dx + dy * dy <= 36) {
                        u(i, j) = std::min(u(i, j), 0.14);
                    }
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    Params p;
    parse_args(argc, argv, p);

    Field u(p.nx, p.ny, 0.0), v(p.nx, p.ny, 0.0), w(p.nx, p.ny, 0.0);
    Field u_next(p.nx, p.ny, 0.0), v_next(p.nx, p.ny, 0.0), w_next(p.nx, p.ny, 0.0);
    Field visual(p.nx, p.ny, 0.0);
    Field hetero(p.nx, p.ny, 0.0), phase_map(p.nx, p.ny, 0.0);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> noise(-4e-4, 4e-4);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);

    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            u(i, j) = 0.02 + noise(rng);
            v(i, j) = 0.01 + noise(rng);
            w(i, j) = 0.02 + noise(rng);
            if (!in_dish(p, i, j)) {
                u(i, j) = 0.0;
                v(i, j) = 0.0;
                w(i, j) = 0.0;
            }

            const double cx = 0.5 * (p.nx - 1);
            const double cy = 0.5 * (p.ny - 1);
            const double dx = (i - cx) / std::max(1.0, 0.5 * p.nx);
            const double dy = (j - cy) / std::max(1.0, 0.5 * p.ny);
            const double r = std::sqrt(dx * dx + dy * dy);
            const double th = std::atan2(dy, dx);
            const double radial_band = std::cos(6.5 * r + 0.7 * th);
            const double azimuth_band = std::sin(3.0 * th - 4.0 * r);
            hetero(i, j) = p.heterogeneity * (0.55 * radial_band + 0.45 * azimuth_band);
            phase_map(i, j) = th;
        }
    }

    apply_layered_source(u, w, p.nx / 2, p.ny / 2, std::max(10, p.nx / 26), p.layered_rings, p.ring_spacing, 0.9, 0.75);
    if (p.wave_mode == Params::WaveMode::Spiral) {
        seed_spiral_break(p, u, w);
    }

    struct CenterSite {
        int x;
        int y;
        int period;
        int phase;
        int radius;
    };
    std::vector<CenterSite> secondary_sites;
    secondary_sites.reserve(std::max(0, p.secondary_centers));
    for (int trial = 0; trial < p.secondary_centers * 40 && static_cast<int>(secondary_sites.size()) < p.secondary_centers;
         ++trial) {
        const int i = static_cast<int>(uni01(rng) * (p.nx - 1));
        const int j = static_cast<int>(uni01(rng) * (p.ny - 1));
        if (!in_dish(p, i, j)) {
            continue;
        }
        if (std::hypot(i - 0.5 * (p.nx - 1), j - 0.5 * (p.ny - 1)) < 0.12 * p.nx) {
            continue;
        }
        const int period = 140 + static_cast<int>(90.0 * uni01(rng));
        const int phase = static_cast<int>(period * uni01(rng));
        const int radius = 4 + static_cast<int>(7.0 * uni01(rng));
        secondary_sites.push_back({i, j, period, phase, radius});
    }

    std::unique_ptr<X11Renderer> renderer;
    if (p.display_every > 0) {
        try {
            renderer = std::make_unique<X11Renderer>(p);
        } catch (const std::exception& e) {
            std::cerr << "[警告] " << e.what() << "\n";
            std::cerr << "自动切换为无界面模式（可加 --save 保存帧）\n";
            p.display_every = 0;
        }
    }

    int frame_id = 0;
    if (p.save_frames) {
        build_visual_field(p, u, w, visual);
        write_ppm(p, visual, frame_id++, p.output_dir);
    }

    auto last_draw = std::chrono::steady_clock::now();

    for (int step = 1; step <= p.steps; ++step) {
        for (int j = 0; j < p.ny; ++j) {
            for (int i = 0; i < p.nx; ++i) {
                if (!in_dish(p, i, j)) {
                    u_next(i, j) = 0.0;
                    v_next(i, j) = 0.0;
                    w_next(i, j) = 0.0;
                    continue;
                }

                const double uu = u(i, j);
                const double vv = v(i, j);
                const double ww = w(i, j);

                const double theta = phase_map(i, j);
                const double eps_local = std::max(0.01, p.epsilon * (1.0 + 0.25 * hetero(i, j)));
                const double phi_local = std::max(
                    0.01, p.phi * (1.0 + 0.35 * hetero(i, j) + p.temporal_mod * std::sin(0.0012 * step + 3.0 * theta)));

                const double reaction_u = (uu * (1.0 - uu) - p.f * vv * ((uu - p.q) / (uu + p.q))) / eps_local;
                const double reaction_v = uu - vv;
                const double reaction_w = phi_local * (uu - ww);

                const double du = reaction_u + p.Du * laplacian(u, i, j, p.dx);
                const double dv = reaction_v + p.Dv * laplacian(v, i, j, p.dx);
                const double dw = reaction_w + p.Dw * laplacian(w, i, j, p.dx);

                u_next(i, j) = std::clamp(uu + p.dt * du, 0.0, 1.5);
                v_next(i, j) = std::clamp(vv + p.dt * dv, 0.0, 1.5);
                w_next(i, j) = std::clamp(ww + p.dt * dw, 0.0, 1.5);
            }
        }

        std::swap(u.data, u_next.data);
        std::swap(v.data, v_next.data);
        std::swap(w.data, w_next.data);

        if (p.pacemaker && step % p.pacemaker_period == 0) {
            apply_pacemaker(p, step, u, w);
        }
        if (p.secondary_centers > 0 && step >= p.secondary_start) {
            for (const auto& s : secondary_sites) {
                if (((step + s.phase) % s.period) == 0) {
                    const double amp = p.secondary_strength * (0.9 + 0.3 * std::sin(0.001 * step + 0.2 * s.x));
                    apply_layered_source(u, w, s.x, s.y, s.radius, 2, p.ring_spacing, amp, 0.85 * amp);
                }
            }
        }

        if (renderer && step % p.display_every == 0) {
            if (!renderer->process_events()) {
                std::cout << "窗口关闭，模拟提前结束。\n";
                break;
            }

            build_visual_field(p, u, w, visual);
            renderer->draw_visual(visual, step);

            auto now = std::chrono::steady_clock::now();
            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_draw).count();
            if (elapsed_ms < 16) {
                std::this_thread::sleep_for(std::chrono::milliseconds(16 - elapsed_ms));
            }
            last_draw = std::chrono::steady_clock::now();
        }

        if (p.save_frames && step % p.frame_every == 0) {
            build_visual_field(p, u, w, visual);
            write_ppm(p, visual, frame_id++, p.output_dir);
        }

        if (step % 1000 == 0) {
            std::cout << "Completed step " << step << " / " << p.steps << "\n";
        }
    }

    std::cout << "Done.\n";
    if (p.save_frames) {
        std::cout << "Frames written to ./'" << p.output_dir << "'\n";
        std::cout << "ffmpeg -framerate 30 -i " << p.output_dir
                  << "/frame_%06d.ppm -pix_fmt yuv420p bz.mp4\n";
    }

    return 0;
}

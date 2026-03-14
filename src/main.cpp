#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <random>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

struct Params {
    int nx = 320;
    int ny = 320;
    double dx = 1.0;
    double dt = 0.0022;
    int steps = 260000;
    double dt = 0.0032;
    int steps = 260000;
    int nx = 260;
    int ny = 260;
    double dx = 1.0;
    double dt = 0.0035;
    int steps = 240000;
    int display_every = 2;

    double epsilon = 0.04;
    double f = 1.35;
    double q = 0.002;
    double phi = 0.075;

    double Du = 0.55;
    double Dv = 0.32;
    double Dw = 0.10;
    double Du = 1.0;
    double Dv = 0.55;
    double Dw = 0.15;
    int nx = 240;
    int ny = 240;
    double dx = 1.0;
    double dt = 0.004;
    int steps = 200000;
    int display_every = 2;

    double epsilon = 0.04;
    double f = 1.4;
    double q = 0.002;
    double phi = 0.08;

    double Du = 1.0;
    double Dv = 0.6;
    double Dw = 0.2;

    bool save_frames = false;
    int frame_every = 20;
    std::string output_dir = "output";
    int pixel_scale = 2;

    bool pacemaker = true;
    int pacemaker_period = 260;
    int pacemaker_radius = 10;
    int layered_rings = 3;
    int ring_spacing = 12;
    int pacemaker_period = 150;
    int pixel_scale = 3;

    bool pacemaker = true;
    int pacemaker_period = 180;
    int pacemaker_radius = 7;

    enum class ViewField { U, W, Mix };
    ViewField view_field = ViewField::Mix;

    enum class WaveMode { Outward, Inward, Dual, Spiral };
    WaveMode wave_mode = WaveMode::Spiral;

    bool auto_contrast = true;
    bool dish_render = true;
    enum class WaveMode { Outward, Inward, Dual };
    WaveMode wave_mode = WaveMode::Outward;

    bool auto_contrast = true;
    int pacemaker_period = 320;
    int pacemaker_radius = 6;

    enum class ViewField { U, W, Mix };
    ViewField view_field = ViewField::Mix;
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

static void map_ferroin_color(double value, unsigned char& r, unsigned char& g, unsigned char& b) {
    double x = std::clamp(value, 0.0, 1.0);
    r = static_cast<unsigned char>(255.0 * (1.0 - x));
    g = static_cast<unsigned char>(70.0 * (1.0 - std::abs(2.0 * x - 1.0)));
    g = static_cast<unsigned char>(60.0 * (1.0 - std::abs(2.0 * x - 1.0)));
    b = static_cast<unsigned char>(255.0 * x);
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
            double v = view_value(p, u(i, j), w(i, j));
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
    for (double& x : visual.data) {
        x = std::clamp((x - vmin) * inv, 0.0, 1.0);
    }
}

static void write_ppm(const Field& visual, int frame_id, const std::string& out_dir) {
    // Mix gives clearer moving excitation fronts for beginners.
    return std::clamp(0.7 * u + 0.3 * w, 0.0, 1.0);
}

static void write_ppm(const Field& w, int frame_id, const std::string& out_dir) {
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
            unsigned char r, g, b;
            map_ferroin_color(visual(i, j), r, g, b);
    ofs << "P6\n" << w.nx << " " << w.ny << "\n255\n";

    for (int j = 0; j < w.ny; ++j) {
        for (int i = 0; i < w.nx; ++i) {
            unsigned char r, g, b;
            map_ferroin_color(w(i, j), r, g, b);
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

        XStoreName(display_, window_, "BZ Reaction: red (reduced) -> blue (oxidized)");
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
            if (ev.type == DestroyNotify) {
                return false;
            }
            if (ev.type == KeyPress) {
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
                map_ferroin_color(visual(i, j), r, g, b);
    void draw_field(const Field& u, const Field& w, int step) {
        for (int j = 0; j < params_.ny; ++j) {
            for (int i = 0; i < params_.nx; ++i) {
                unsigned char r, g, b;
                map_ferroin_color(view_value(params_, u(i, j), w(i, j)), r, g, b);
                unsigned long pixel = channel_to_masked(r, red_mask_) | channel_to_masked(g, green_mask_) |
                                      channel_to_masked(b, blue_mask_);

                for (int sy = 0; sy < params_.pixel_scale; ++sy) {
                    const int y = j * params_.pixel_scale + sy;
                    auto* row = reinterpret_cast<std::uint32_t*>(image_buffer_.data() + y * image_stride_);
                    for (int sx = 0; sx < params_.pixel_scale; ++sx) {
                        const int x = i * params_.pixel_scale + sx;
                    int y = j * params_.pixel_scale + sy;
                    auto* row = reinterpret_cast<std::uint32_t*>(image_buffer_.data() + y * image_stride_);
                    for (int sx = 0; sx < params_.pixel_scale; ++sx) {
                        int x = i * params_.pixel_scale + sx;
                        row[x] = static_cast<std::uint32_t>(pixel);
                    }
                }
            }
        }

        XPutImage(display_, window_, gc_, image_, 0, 0, 0, 0, width_, height_);
        std::ostringstream title;
        title << "BZ spiral wave step=" << step << " (按任意键退出)";
        title << "BZ Reaction step=" << step << " (按任意键退出)";
        title << "BZ Reaction  step=" << step << " (按任意键退出)";
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
            p.nx = std::max(16, std::stoi(a.substr(5)));
        } else if (a.rfind("--ny=", 0) == 0) {
            p.ny = std::max(16, std::stoi(a.substr(5)));
        } else if (a == "--no-pacemaker") {
            p.pacemaker = false;
        } else if (a.rfind("--pacemaker-period=", 0) == 0) {
            p.pacemaker_period = std::max(20, std::stoi(a.substr(19)));
        } else if (a.rfind("--layers=", 0) == 0) {
            p.layered_rings = std::max(1, std::stoi(a.substr(9)));
        } else if (a.rfind("--ring-spacing=", 0) == 0) {
            p.ring_spacing = std::max(4, std::stoi(a.substr(15)));
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
            const std::string v = a.substr(11);
            p.auto_contrast = (v != "fixed");
        }
    }
}

static void apply_source_patch(Field& u, Field& w, int cx, int cy, int radius, double u_peak, double w_peak) {
    for (int j = 0; j < u.ny; ++j) {
        for (int i = 0; i < u.nx; ++i) {
            const int dx = i - cx;
            const int dy = j - cy;
            int dx = i - cx;
            int dy = j - cy;
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

static void apply_pacemaker(const Params& p, Field& u, Field& w) {
    if (!p.pacemaker) {
        return;
    }

    const int cx = p.nx / 2;
    const int cy = p.ny / 2;

    if (p.wave_mode == Params::WaveMode::Outward || p.wave_mode == Params::WaveMode::Dual) {
        apply_layered_source(u, w, cx, cy, p.pacemaker_radius, p.layered_rings, p.ring_spacing, 0.85, 0.72);
        apply_source_patch(u, w, cx, cy, p.pacemaker_radius, 1.1, 0.9);
    }

    if (p.wave_mode == Params::WaveMode::Inward || p.wave_mode == Params::WaveMode::Dual) {
        const int ex = p.nx / 2;
        const int ey = p.ny / 10;
        apply_source_patch(u, w, ex, ey, p.pacemaker_radius, 0.95, 0.78);
        apply_source_patch(u, w, p.nx - ex, ey, p.pacemaker_radius, 0.95, 0.78);
        apply_source_patch(u, w, ex, p.ny - ey, p.pacemaker_radius, 0.95, 0.78);
        apply_source_patch(u, w, p.nx - ex, p.ny - ey, p.pacemaker_radius, 0.95, 0.78);
    }

    if (p.wave_mode == Params::WaveMode::Spiral) {
        apply_layered_source(u, w, p.nx / 3, p.ny / 2, p.pacemaker_radius, p.layered_rings, p.ring_spacing, 0.82,
                             0.70);
        apply_layered_source(u, w, 2 * p.nx / 3, p.ny / 2, p.pacemaker_radius, p.layered_rings, p.ring_spacing,
                             0.82, 0.70);
        const int ey = p.ny / 12;
        apply_source_patch(u, w, ex, ey, p.pacemaker_radius, 1.05, 0.85);
        apply_source_patch(u, w, p.nx - ex, ey, p.pacemaker_radius, 1.05, 0.85);
        apply_source_patch(u, w, ex, p.ny - ey, p.pacemaker_radius, 1.05, 0.85);
        apply_source_patch(u, w, p.nx - ex, p.ny - ey, p.pacemaker_radius, 1.05, 0.85);
    }

    if (p.wave_mode == Params::WaveMode::Spiral) {
        apply_source_patch(u, w, p.nx / 3, p.ny / 2, p.pacemaker_radius, 1.07, 0.88);
        apply_source_patch(u, w, 2 * p.nx / 3, p.ny / 2, p.pacemaker_radius, 1.07, 0.88);
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

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> noise(-4e-4, 4e-4);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> noise(-2e-4, 2e-4);

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
        }
    }

    apply_layered_source(u, w, p.nx / 2, p.ny / 2, std::max(10, p.nx / 26), p.layered_rings, p.ring_spacing, 0.9, 0.75);
    if (p.wave_mode == Params::WaveMode::Spiral) {
        seed_spiral_break(p, u, w);
    apply_source_patch(u, w, p.nx / 2, p.ny / 2, std::max(8, p.nx / 28), 1.0, 0.85);
    if (p.wave_mode == Params::WaveMode::Spiral) {
        seed_spiral_break(p, u, w);
        }
    }

    // Initial trigger
    apply_source_patch(u, w, p.nx / 2, p.ny / 2, std::max(8, p.nx / 28), 1.0, 0.85);
    int cx = p.nx / 3;
    int cy = p.ny / 2;
    int radius = std::max(8, p.nx / 25);
    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            int dx = i - cx;
            int dy = j - cy;
            if (dx * dx + dy * dy <= radius * radius) {
                u(i, j) = 0.85;
                v(i, j) = 0.20;
                w(i, j) = 0.80;
            }
        }
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
        Field visual(p.nx, p.ny, 0.0);
        for (int j = 0; j < p.ny; ++j) {
            for (int i = 0; i < p.nx; ++i) {
                visual(i, j) = view_value(p, u(i, j), w(i, j));
            }
        }
        write_ppm(visual, frame_id++, p.output_dir);
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

                const double reaction_u = (uu * (1.0 - uu) - p.f * vv * ((uu - p.q) / (uu + p.q))) / p.epsilon;
                const double reaction_v = uu - vv;
                const double reaction_w = p.phi * (uu - ww);

                const double du = reaction_u + p.Du * laplacian(u, i, j, p.dx);
                const double dv = reaction_v + p.Dv * laplacian(v, i, j, p.dx);
                const double dw = reaction_w + p.Dw * laplacian(w, i, j, p.dx);
                double uu = u(i, j);
                double vv = v(i, j);
                double ww = w(i, j);

                double reaction_u = (uu * (1.0 - uu) - p.f * vv * ((uu - p.q) / (uu + p.q))) / p.epsilon;
                double reaction_v = uu - vv;
                double reaction_w = p.phi * (uu - ww);

                double du = reaction_u + p.Du * laplacian(u, i, j, p.dx);
                double dv = reaction_v + p.Dv * laplacian(v, i, j, p.dx);
                double dw = reaction_w + p.Dw * laplacian(w, i, j, p.dx);

                u_next(i, j) = std::clamp(uu + p.dt * du, 0.0, 1.5);
                v_next(i, j) = std::clamp(vv + p.dt * dv, 0.0, 1.5);
                w_next(i, j) = std::clamp(ww + p.dt * dw, 0.0, 1.5);
            }
        }

        std::swap(u.data, u_next.data);
        std::swap(v.data, v_next.data);
        std::swap(w.data, w_next.data);

        if (p.pacemaker && step % p.pacemaker_period == 0) {
            apply_pacemaker(p, u, w);
            const int px = p.nx / 4;
            const int py = p.ny / 2;
            for (int j = 0; j < p.ny; ++j) {
                for (int i = 0; i < p.nx; ++i) {
                    const int dx = i - px;
                    const int dy = j - py;
                    if (dx * dx + dy * dy <= p.pacemaker_radius * p.pacemaker_radius) {
                        u(i, j) = 1.1;
                        w(i, j) = std::max(w(i, j), 0.9);
                    }
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
            renderer->draw_field(u, w, step);

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
            Field visual(p.nx, p.ny, 0.0);
            for (int j = 0; j < p.ny; ++j) {
                for (int i = 0; i < p.nx; ++i) {
                    visual(i, j) = view_value(p, u(i, j), w(i, j));
                }
            }
            write_ppm(visual, frame_id++, p.output_dir);
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

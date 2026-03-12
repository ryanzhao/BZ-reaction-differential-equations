#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

struct Params {
    int nx = 220;
    int ny = 220;
    double dx = 1.0;
    double dt = 0.005;
    int steps = 12000;
    int frame_every = 40;

    double epsilon = 0.04;
    double f = 1.4;
    double q = 0.002;
    double phi = 0.08;

    double Du = 1.0;
    double Dv = 0.6;
    double Dw = 0.2;
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
    // Neumann-like boundary handling by mirrored edge indexing.
    int il = clamp_idx(i - 1, 0, a.nx - 1);
    int ir = clamp_idx(i + 1, 0, a.nx - 1);
    int jb = clamp_idx(j - 1, 0, a.ny - 1);
    int jt = clamp_idx(j + 1, 0, a.ny - 1);

    return (a(il, j) + a(ir, j) + a(i, jb) + a(i, jt) - 4.0 * a(i, j)) / (dx * dx);
}

static void write_ppm(const Field& w, int frame_id, const std::string& out_dir) {
    std::filesystem::create_directories(out_dir);

    std::ostringstream name;
    name << out_dir << "/frame_" << std::setw(5) << std::setfill('0') << frame_id << ".ppm";

    std::ofstream ofs(name.str(), std::ios::binary);
    ofs << "P6\n" << w.nx << " " << w.ny << "\n255\n";

    for (int j = 0; j < w.ny; ++j) {
        for (int i = 0; i < w.nx; ++i) {
            double x = std::clamp(w(i, j), 0.0, 1.0);
            // Ferroin-inspired palette: reduced(red) -> oxidized(blue).
            unsigned char r = static_cast<unsigned char>(255.0 * (1.0 - x));
            unsigned char g = static_cast<unsigned char>(60.0 * (1.0 - std::abs(2.0 * x - 1.0)));
            unsigned char b = static_cast<unsigned char>(255.0 * x);
            ofs.write(reinterpret_cast<char*>(&r), 1);
            ofs.write(reinterpret_cast<char*>(&g), 1);
            ofs.write(reinterpret_cast<char*>(&b), 1);
        }
    }
}

int main() {
    Params p;

    Field u(p.nx, p.ny, 0.0);
    Field v(p.nx, p.ny, 0.0);
    Field w(p.nx, p.ny, 0.0);

    Field u_next(p.nx, p.ny, 0.0);
    Field v_next(p.nx, p.ny, 0.0);
    Field w_next(p.nx, p.ny, 0.0);

    // Rest state + low-amplitude noise.
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> noise(-2e-4, 2e-4);

    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            u(i, j) = 0.02 + noise(rng);
            v(i, j) = 0.01 + noise(rng);
            w(i, j) = 0.02 + noise(rng);
        }
    }

    // Seed a circular perturbation for a target wave.
    int cx = p.nx / 3;
    int cy = p.ny / 2;
    int radius = 10;
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

    int frame_id = 0;
    write_ppm(w, frame_id++, "output");

    for (int step = 1; step <= p.steps; ++step) {
        for (int j = 0; j < p.ny; ++j) {
            for (int i = 0; i < p.nx; ++i) {
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

        if (step % p.frame_every == 0) {
            write_ppm(w, frame_id++, "output");
        }

        if (step % 1000 == 0) {
            std::cout << "Completed step " << step << " / " << p.steps << "\n";
        }
    }

    std::cout << "Done. Frames written to ./output\n";
    std::cout << "Example video command: ffmpeg -framerate 25 -i output/frame_%05d.ppm -pix_fmt yuv420p bz.mp4\n";

    return 0;
}

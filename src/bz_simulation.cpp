#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct Params {
    int nx = 200;
    int ny = 200;
    int steps = 8000;
    int output_every = 100;

    double dx = 1.0;
    double dt = 0.001;

    // Oregonator-like parameters (dimensionless)
    double epsilon = 0.05;
    double q = 0.002;
    double f = 1.2;
    double delta = 0.02;

    // Diffusion coefficients
    double Dx = 1.0;
    double Dy = 0.6;
    double Dz = 0.0;
};

struct Fields {
    std::vector<double> X, Y, Z;
    std::vector<double> Xnext, Ynext, Znext;
};

int idx(int i, int j, int nx) { return j * nx + i; }

int clamp_mirror(int a, int n) {
    if (a < 0) return 1;         // mirror for Neumann-like boundary
    if (a >= n) return n - 2;
    return a;
}

double laplacian(const std::vector<double>& u, int i, int j, int nx, int ny, double dx) {
    const int ip = clamp_mirror(i + 1, nx);
    const int im = clamp_mirror(i - 1, nx);
    const int jp = clamp_mirror(j + 1, ny);
    const int jm = clamp_mirror(j - 1, ny);

    return (u[idx(ip, j, nx)] + u[idx(im, j, nx)] + u[idx(i, jp, nx)] + u[idx(i, jm, nx)] -
            4.0 * u[idx(i, j, nx)]) /
           (dx * dx);
}

void initialize(const Params& p, Fields& flds) {
    const int n = p.nx * p.ny;
    flds.X.assign(n, 0.1);
    flds.Y.assign(n, 0.1);
    flds.Z.assign(n, 0.1);
    flds.Xnext.resize(n);
    flds.Ynext.resize(n);
    flds.Znext.resize(n);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> noise(-0.005, 0.005);

    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            int k = idx(i, j, p.nx);
            flds.X[k] += noise(rng);
            flds.Y[k] += noise(rng);
            flds.Z[k] += noise(rng);
        }
    }

    // Trigger pulse near center to initiate target waves.
    const int cx = p.nx / 2;
    const int cy = p.ny / 2;
    const int r2 = 10 * 10;
    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            int dx = i - cx;
            int dy = j - cy;
            if (dx * dx + dy * dy < r2) {
                int k = idx(i, j, p.nx);
                flds.X[k] = 0.8;
                flds.Z[k] = 0.8;
            }
        }
    }
}

void write_ppm(const Params& p, const Fields& flds, int frame_id) {
    std::ostringstream name;
    name << "frame_" << std::setfill('0') << std::setw(5) << frame_id << ".ppm";
    std::ofstream out(name.str(), std::ios::binary);
    out << "P6\n" << p.nx << " " << p.ny << "\n255\n";

    // Color map based on Z (oxidized fraction proxy): red -> blue
    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            const double z = std::clamp(flds.Z[idx(i, j, p.nx)], 0.0, 1.0);
            const uint8_t r = static_cast<uint8_t>((1.0 - z) * 220 + 20);
            const uint8_t g = static_cast<uint8_t>((1.0 - std::abs(z - 0.5) * 2.0) * 40 + 10);
            const uint8_t b = static_cast<uint8_t>(z * 230 + 20);
            out.write(reinterpret_cast<const char*>(&r), 1);
            out.write(reinterpret_cast<const char*>(&g), 1);
            out.write(reinterpret_cast<const char*>(&b), 1);
        }
    }
}

void step(const Params& p, Fields& flds) {
    for (int j = 0; j < p.ny; ++j) {
        for (int i = 0; i < p.nx; ++i) {
            const int k = idx(i, j, p.nx);
            const double X = flds.X[k];
            const double Y = flds.Y[k];
            const double Z = flds.Z[k];

            const double rx = (p.q * Y - X * Y + X * (1.0 - X)) / p.epsilon;
            const double ry = -p.q * Y - X * Y + p.f * Z;
            const double rz = p.delta * (X - Z);

            const double diff_x = p.Dx * laplacian(flds.X, i, j, p.nx, p.ny, p.dx);
            const double diff_y = p.Dy * laplacian(flds.Y, i, j, p.nx, p.ny, p.dx);
            const double diff_z = p.Dz * laplacian(flds.Z, i, j, p.nx, p.ny, p.dx);

            flds.Xnext[k] = X + p.dt * (rx + diff_x);
            flds.Ynext[k] = Y + p.dt * (ry + diff_y);
            flds.Znext[k] = Z + p.dt * (rz + diff_z);

            // Keep fields bounded for numerical robustness.
            flds.Xnext[k] = std::clamp(flds.Xnext[k], 0.0, 2.0);
            flds.Ynext[k] = std::clamp(flds.Ynext[k], 0.0, 2.0);
            flds.Znext[k] = std::clamp(flds.Znext[k], 0.0, 1.0);
        }
    }

    std::swap(flds.X, flds.Xnext);
    std::swap(flds.Y, flds.Ynext);
    std::swap(flds.Z, flds.Znext);
}

int main() {
    Params p;
    Fields flds;
    initialize(p, flds);

    int frame_id = 0;
    write_ppm(p, flds, frame_id++);

    for (int n = 1; n <= p.steps; ++n) {
        step(p, flds);
        if (n % p.output_every == 0) {
            write_ppm(p, flds, frame_id++);
            std::cout << "step " << n << "/" << p.steps << "\n";
        }
    }

    std::cout << "Done. Convert frames to video, e.g.:\n"
              << "ffmpeg -framerate 20 -i frame_%05d.ppm -pix_fmt yuv420p bz.mp4\n";
    return 0;
}

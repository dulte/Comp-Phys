#include "harmonicoscillator.h"

HarmonicOscillator::HarmonicOscillator(double omega_) {
    omega = omega_;
}

double HarmonicOscillator::computePotential(double r) {
    return omega*omega*r*r;
}

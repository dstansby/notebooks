"""
Based on
  http://www.physics.sfasu.edu/astro/color/spectra.html
Taken from
  https://hsugawa.blogspot.co.uk/2010/01/matplotlib-colormap-for-visible.html

  RGB VALUES FOR VISIBLE WAVELENGTHS   by Dan Bruton (astro@tamu.edu)
"""

import numpy as np


def factor(wl):
    return np.select(
        [ wl > 700.,
          wl < 420.,
          True],
        [ .3 + .7 * (780. - wl) / (780. - 700.),
          .3 + .7 * (wl - 380.) / (420. - 380.),
          1.0])


def raw_r(wl):
    return np.select(
        [ wl >= 580.,
          wl >= 510.,
          wl >= 440.,
          wl >= 380.,
          True],
        [ 1.0,
          (wl - 510.) / (580. - 510.),
          0.0,
          (wl - 440.) / (380. - 440.),
          0.0])


def raw_g(wl):
    return np.select(
        [ wl >= 645.,
          wl >= 580.,
          wl >= 490.,
          wl >= 440.,
          True ],
        [ 0.0,
          (wl - 645.) / (580. - 645.),
          1.0,
          (wl - 440.) / (490. - 440.),
          0.0])


def raw_b(wl):
    return np.select(
        [ wl >= 510.,
          wl >= 490.,
          wl >= 380.,
          True],
        [ 0.0,
          (wl - 510.) / (490. - 510.),
          1.0,
          0.0])


gamma = 0.80


def correct_r(wl):
    return np.power(factor(wl) * raw_r(wl), gamma)


def correct_g(wl):
    return np.power(factor(wl) * raw_g(wl), gamma)


def correct_b(wl):
    return np.power(factor(wl) * raw_b(wl), gamma)


def rgb(wl):
    '''
    Given wavelength in nano-metres, return RGB values.
    '''
    return np.transpose([correct_r(wl), correct_g(wl), correct_b(wl)])

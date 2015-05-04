# Introduction #

This is a Python package for computing the electromagnetic scattering properties of nonspherical particles using the _T_-matrix method.

The package can be used particles of any size as long as the wavelength is not too long compared to the particle size. However, there is some functionality that is intended specifically for the calculation of scattering properties of hydrometeors at microwave frequencies, such as those used in weather radars.



# Installation #
Make sure you have [NumPy](http://numpy.org/) and [SciPy](http://scipy.org/) (and of course, a Python environment) installed. You also need a Fortran 77 compiler (e.g. gfortran) to install the package. For example, if you are on a system with `apt`, this should get you what you need:
```
$ sudo apt-get install gfortran python-numpy python-scipy
```

The code is easiest to install from PyPI. If you have  [pip](http://pypi.python.org/pypi/pip) installed, install simply with
```
$ pip install pytmatrix
```
(using sudo if installing globally).

If you can't use `pip`, install manually:
  1. [Download](http://code.google.com/p/pytmatrix/downloads/list) the latest package and unzip OR [clone](http://code.google.com/p/pytmatrix/source/checkout) the latest version from the repository.
  1. Go to the root directory of the package and run (with sudo if you need the permissions):
```
$ python setup.py install 
```

# Testing #

This is to verify that you have a working package after installation. Open a Python shell. Run the following:
```
from pytmatrix.test import test_tmatrix
test_tmatrix.run_tests()
```
An output like this should be printed:
```
test_adaptive_orient (pytmatrix.test.test_tmatrix.TMatrixTests)
Test an adaptive orientation averaging case ... ok
test_asymmetry (pytmatrix.test.test_tmatrix.TMatrixTests)
Test calculation of the asymmetry parameter ... ok
test_fixed_orient (pytmatrix.test.test_tmatrix.TMatrixTests)
Test a fixed-point orientation averaging case ... ok
test_optical_theorem (pytmatrix.test.test_tmatrix.TMatrixTests)
Optical theorem: test that for a lossless particle, Csca=Cext ... ok
test_psd (pytmatrix.test.test_tmatrix.TMatrixTests)
Test a case that integrates over a particle size distribution ... ok
test_radar (pytmatrix.test.test_tmatrix.TMatrixTests)
Test that the radar properties are computed correctly ... ok
test_rayleigh (pytmatrix.test.test_tmatrix.TMatrixTests)
Test match with Rayleigh scattering for small spheres ... ok
test_single (pytmatrix.test.test_tmatrix.TMatrixTests)
Test a single-orientation case ... ok

----------------------------------------------------------------------
Ran 8 tests in 12.072s

OK
```
The running time may vary; the most important thing is that "ok" is printed after each test. If some of the tests fail, please notify the author.

Please note that these tests only verify the basic functionality of the package, they do not (as of the latest version) test for convergence in edge cases.

# Usage #

## Minimal example ##
This calculates the amplitude scattering matrix **S** for a sphere with a volume-equivalent radius of 1.0 and horizontal-to-vertical axis ratio of 0.6 and refractive index 1.5+0.5i, at a wavelength of 10.0:
```
>>> from pytmatrix import tmatrix
>>> scatterer = tmatrix.Scatterer(radius=1.0, wavelength=10.0, m=complex(1.5, 0.5), axis_ratio=1.0/0.6)
>>> scatterer.get_S()
array([[  9.66435719e-02 +6.79614074e-02j,
          6.16862803e-25 +7.07135826e-25j],
       [ -6.16862447e-25 -7.07128493e-25j,
         -1.01600111e-01 -1.06748868e-01j]])
```

## Using the `Scatterer` class ##

The `Scatterer` class, found in the `pytmatrix.tmatrix` module, is an object-oriented interface to the Fortran 77 _T_-matrix code. To use it, create an instance of the class, set the properties of the scatterer, and then run one of the commands to retrieve the amplitude and/or phase matrices.

### Attributes ###

The properties of the scatterer, the radiation and the computations can be set using the attributes of the `TMatrix` instance. You can also set any of them as keyword arguments to the constructor, as shown in the example above. The available properties are:
| **Attribute** | **Description** | **Default value** |
|:--------------|:----------------|:------------------|
| `radius` | The radius of the equivalent sphere (see _rat_). | `1.0` |
| `rat` | If this is set to 1, `radius` is the equivalent-volume radius, otherwise `radius` is the equivalent-surface-area radius. | `1.0` |
| `wavelength` | The wavelength of the incident radiation. | `1.0` |
| `m` | The complex refractive index of the scatterer. | `complex(2,0)` |
| `axis_ratio` | The horizontal-to-rotational axis ratio. | `1.000001` |
| `np` | The type of the particles. Use `np=-1` for spheroids and `np=-2` for cylinders. Note that Chebyshev particles are not currently supported. | `-1` |
| `ddelt` | The accuracy of the computations. | `1e-3` |
| `ndgs` | Number of division points used to integrate over the particle surface. Increase this if the computations do not converge. | `2` |
| `alpha` | The Euler angle α (alpha) of the scatterer orientation. | `0.0` |
| `beta` | The Euler angle β (beta) of the scatterer orientation. | `0.0` |
| `thet0` | The zenith angle of the incident beam. | `90.0` |
| `thet` | The zenith angle of the scattered beam. | `90.0` |
| `phi0` | The azimuth angle of the incident beam. | `0.0` |
| `phi` | The azimuth angle of the scattered beam. | `180.0` |
| `Kw_sqr` | The reference value of the water dielectric factor _|K|_² (used for calculating radar reflectivity). | `0.93` |
| `orient` | The method of performing orientation averaging. | `orientation.orient_single` |
| `or_pdf` | Probability distribution function for the particle orientation. | `orientation.gaussian_pdf()` |
| `n_alpha` | Number of integration points for averaging over the _alpha_ angle when fixed-point orientation averaging is used. | `5` |
| `n_beta` | As above, but for the `beta` angle. | `10` |
| `psd_integrator` | Set this to a PSDIntegrator instance to enable size distribution integration. If this is `None`, size distribution integration is not used. See the PSD integration documentation for more information. | `None` |
| `psd` | Set to a callable object giving the PSD value for a given diameter (for example a GammaPSD instance); default None. Has no effect if psd\_integrator is `None`. | `None` |

**Note on units**: `radius` and `wavelength` should always have the same unit. Some components of the `radar` module expect them to be given in millimeters, but otherwise, other units can be used. The results of the scattering computations will be in corresponding units. Angles should always be given in _degrees_ (not radians).

The names of the parameters are defined such as to be compatible with the original _T_-Matrix code.

### Methods ###

The `Scatterer` class encapsulates the computation of the amplitude (**S**) and phase (**Z**) matrices of the scatterer. There are three instance methods that can be used to access the results, and an additional convenience function for setting the geometry.

| **Method** | **Arguments** | **Description** |
|:-----------|:--------------|:----------------|
| `get_SZ()` | None | Returns the amplitude and phase matrices for the currently defined scatterer, packed in a tuple. |
| `get_S()` | None | As get\_SZ, but only returns the amplitude matrix. |
| `get_Z()` | None | As get\_SZ, but only returns the phase matrix. |
| `set_geometry(geom)` | `geom`: A tuple of (`thet0`, `thet`, `phi0`, `phi`, `alpha`, `beta`). | A convenience function for setting all angles of the scattering geometry at the same time (see attributes). |

One of the advantages of the _T_-matrix method is that once the _T_ matrix is computed, it is relatively inexpensive to compute the amplitude and phase matrices for any scattering geometry. The `Scatterer` class implements a lazy evaluation scheme that exploits this property. Calling any of `get_SZ`, `get_S` or `get_Z` will reuse the results of previous calculations if they are still applicable. Changing any of the following parameters will trigger a recalculation of the _T_ matrix : `radius`, `rat`, `wavelength`, `m`, `axis_ratio`, `np`, `ddelt`, `ndgs`. The recalculation is not done immediately, but instead on the next call to one of the above-mentioned methods. Changing any of the other attributes will only recalculate the amplitude and phase matrices. If no parameters are changed between calls to those methods, the amplitude and phase matrices will also be reused.

## Orientation averaging: the `orientation` module ##

This module contains functions that allow the use of different methods for computing the scattering properties. Most importantly, they allow the use of orientation averaging when computing the scattering properties.

To enable orientation averaging, two things are needed. Firstly, you need to set the `orient` property of your `Scatterer` object to a function that is able to perform averaging. The function `orientation.orient_single`, the default, only computes for a single orientation. For orientation averaging, set `orient` to either `orientation.orient_averaged_fixed` or `orientation.orient_averaged_adaptive`. The former uses a weighted fixed-point quadrature integration for averaging, while the latter, much slower, employs an adaptive integration method. Both orientation averaging methods fully utilize the possibility to reuse the _T_ matrix. If you enable orientation averaging, the `Scatterer` attributes _alpha_ and _beta_ are ignored.

Secondly, you should set the orientation distribution function. You can either use a function from the `orientation` module, or define your own. Currently, two different distributions are implemented: a uniform distribution returned by `orientation.uniform_pdf()` and a Gaussian distribution in `orientation.gaussian_pdf()`. The latter takes one argument, which specifies the standard deviation of the angle with respect to vertical orientation (this is often called the canting angle).

If you define your own orientation distribution, it should be represented by a callable object (e.g. a function) which takes one argument, the deviation from vertical orientation in degrees, and returns the orientation probability for that angle. Take a look at the functions in the `orientation` module for an example. The distribution function must be defined in the interval [0,180] and must integrate to 1 over that interval. These criteria are not automatically checked, so you are responsible for the correct behavior of your function.

## Particle size distribution integration: the `psd` module ##

To enable integration over particle size distributions (PSD), you should set the `psd_integrator` and `psd` attributes of the `Scatterer` instance. When PSDs are used, it is common to use the same particles for several different PSDs. For this reason, the `PSDIntegrator` class computes and stores a lookup table of the amplitude and phase matrices over a particle size interval; the PSD-weighted integration is then performed over this table.

The `psd_integrator` attribute should be set to an instance of the `PSDIntegrator` class from the `psd` module. The following attributes can be set:
| **Attribute** | **Description** | **Default value** |
|:--------------|:----------------|:------------------|
| `num_points` | The number of different (equally spaced) particle diameters at which to store the amplitude and phase matrices. | 500 |
| `D_max` | The maximum diameter for which to store the amplitude and phase matrices. | `psd.Dmax` |
| `m_func` | The refractive index as a function of size (see below). | `None` |
| `axis_ratio_func` | The horizontal-to-rotational axis ratio as a function of size (see below). | `None` |
| `geometries` | A tuple of the scattering geometries for which to store the amplitude and phase matrices. The format of a geometry specification is the same as that for the `set_geometry` function. The geometry of the `Scatterer` instance must be set to one of these geometries. | `(tmatrix_aux.geom_horiz_back,)` |

The `m_func` and `axis_ratio_func` attributes specify the refractive index and the axis ratio as a function of size. These should be callables that take the diameter as an argument and return the corresponding variable. They must be defined in the interval (0,`D_max`]. If set to `None`, the constant values of `m` and `axis_ratio` are used, respectively.

The `psd` attribute specifies the particle size distribution. Currently, two PSDs are pre-defined in the `psd` module: the normalized gamma PSD, as given by Eq. (7.62) of Bringi and Chandrasekar (2001), and a binned PSD for computing results for experimental data. The gamma PSD can be constructed as an instance of the `psd.GammaPSD` class, which takes three keyword arguments: `D0`, `Nw` and `mu`. The `psd.BinnedPSD` takes as its arguments the bin edges and bin values.

If you define your own PSD, it should be a callable object with two properties. Firstly, it should be called with a single argument that specifies the diameter, and it should return the value of the PSD for that diameter. The PSD must be defined in the interval (0, `D_max`) and should integrate to the number concentration of particles in a unit volume. Secondly, the PSD may optionally define a way to determine if two PSDs are equal. This avoids reintegrating over the PSD if it has not changed.

To use the `PSDIntegrator` class, set the above attributes as needed, then call the `init_scatter_table` method of the `PSDIntegrator` instance using your `Scatterer` instance as the argument. After this, the PSD-integrated amplitude and phase matrices can be retrieved as usual from the `Scatterer` class.

Because computing the lookup tables can sometimes take a long time, it is possible to save the lookup tables to a file. Call the `save_scatter_table` method to do this and `load_scatter_table` to retrieve the saved tables. After a call to `load_scatter_table`, you can use the `PSDIntegrator` object just as if you had called `init_scatter_table`.

The methods of `PSDIntegrator` instances are summarized below:
| **Method** | **Arguments** | **Description** |
|:-----------|:--------------|:----------------|
| `init_scatter_table()` | None | Initialize the matrix lookup tables |
| `save_scatter_table(fn, desription="")` | `fn`: the name and path of the file to save to<br><code>description</code>: A description the table.  <table><thead><th> Save the matrix lookup tables to a file. </th></thead><tbody>
<tr><td> <code>load_scatter_table(fn)</code> </td><td> <code>fn</code>: the name and path of the file to load from. </td><td> Load the matrix lookup tables from a file. </td></tr></tbody></table>

<h2>The <code>tmatrix_aux</code> module</h2>

This module has some convenient definitions of constants that may make the code easier to use. These are listed below:<br>
<table><thead><th> <b>Constant</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> <code>VERSION</code> </td><td> The current version of <code>pytmatrix</code> </td></tr>
<tr><td> <code>wl_S</code>, <code>wl_C</code>, <code>wl_X</code>, <code>wl_Ka</code>, <code>wl_Ku</code>, <code>wl_Ka</code>, <code>wl_W</code> </td><td> Typical wavelengths (millimeters) for various radar bands. </td></tr>
<tr><td> K_w_sqr </td><td> The reference dielectric factors  <i>|K|</i>² for calculating the radar reflectivity. This is a dictionary with keys for each of the above mentioned wavelengths. </td></tr>
<tr><td> geom_horiz_back </td><td> A geometry tuple in the form accepted by <code>TMatrix.set_geometry</code>. Defines the geometry for backscattering of horizontal incident radiation. </td></tr>
<tr><td> geom_horiz_forw </td><td> As above, but for forward scattering of horizontal incident radiation. </td></tr>
<tr><td> geom_vert_back </td><td> As above, but for backscattering of vertical incident radiation. </td></tr>
<tr><td> geom_vert_forw </td><td> As above, but for forward scattering of vertical incident radiation. </td></tr></tbody></table>

<h2>Computation of scattering parameters: the <code>scatter</code> module</h2>

This module contains functions that are used to compute scattering quantities of general interest in many applications. These all take any <code>Scatterer</code> object as their first argument and compute the results for the selected scatterer parameters. Some functions operate only on specific geometries, which are also described below.<br>
<br>
<table><thead><th> <b>Method</b> </th><th> <b>Arguments</b> </th><th> <b>Description</b> </th><th> <b>Geometry</b> </th></thead><tbody>
<tr><td> <code>sca_intensity(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: a <code>Scatterer</code> object<br><code>h_pol</code>: <code>True</code> for horizontal polarization of incident radiation, <code>False</code> for vertical polarization </td><td> Scattering intensity (phase function) to the direction specified by the <code>thet</code> and <code>phi</code> arguments of <code>scatterer</code>. </td><td> Any </td></tr>
<tr><td> <code>sca_xsect(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: As above. </td><td> Scattering cross section </td><td> Any (<code>thet</code> and <code>phi</code> parameters ignored). </td></tr>
<tr><td> <code>ext_xsect(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: As above. </td><td> Extinction cross section </td><td> Forward scattering (must be set as this is computed using the optical theorem). </td></tr>
<tr><td> <code>ssa(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: As above. </td><td> Single scattering albedo, or <code>sca_xsect/ext_xsect</code>. </td><td> Forward scattering (must be set as the extinction cross section is computed using the optical theorem). </td></tr>
<tr><td> <code>asym(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: As above. </td><td> Asymmetry parameter <i><cos Θ></i> </td><td> Any (<code>thet</code> and <code>phi</code> parameters ignored). </td></tr>
<tr><td> <code>ldr(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: If <code>True</code>, return LDR_h. If <code>False</code>, return LDR_v. </td><td> Linear depolarization ratio </td><td> Backscattering </td></tr></tbody></table>


<h2>Microwave-specific functionality: the <code>radar</code> module</h2>

This module contains functions that are mostly interesting for use with  radars and other microwave instruments. These all take any <code>Scatterer</code>  object as their first argument, although many of them only make sense if the <code>Scatterer</code> object has PSD integration defined by enabling the <code>psd_integrator</code> attribute. Many of the functions also make sense only for a specific scattering geometry; you are responsible for setting an appropriate geometry yourself. Most of these functions expect that you give all sizes in millimeters.<br>
<br>
<table><thead><th> <b>Method</b> </th><th> <b>Arguments</b> </th><th> <b>Description</b> </th><th> <b>Geometry</b> </th></thead><tbody>
<tr><td> <code>radar_xsect(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: a <code>Scatterer</code> object<br><code>h_pol</code>: <code>True</code> for horizontal polarization, <code>False</code> for vertical polarization </td><td> Radar cross section </td><td> Any (use non-backscattering geometries for bistatic cross section) </td></tr>
<tr><td> <code>refl(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: As above </td><td> Radar reflectivity </td><td> Backscattering </td></tr>
<tr><td> <code>ldr(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: If <code>True</code>, return LDR_h. If <code>False</code>, return LDR_v. </td><td> Linear depolarization ratio </td><td> Backscattering </td></tr>
<tr><td> <code>Zdr(scatterer)</code> </td><td> <code>scatterer</code>: As above<br> </td><td> Differential reflectivity </td></tr>
Backscattering <br>
<tr><td> <code>delta_hv(scatterer)</code> </td><td> <code>scatterer</code>: As above<br> </td><td> Scattering differential phase </td><td> Backscattering </td></tr>
<tr><td> <code>rho_hv(scatterer)</code> </td><td> <code>scatterer</code>: As above<br> </td><td> Co-polar correlation </td><td> Backscattering </td></tr>
<tr><td> <code>Kdp(scatterer)</code> </td><td> <code>scatterer</code>: As above<br> </td><td> Specific differential phase (deg/km) </td><td> Forward scattering </td></tr>
<tr><td> <code>Ai(scatterer, h_pol=True)</code> </td><td> <code>scatterer</code>: As above<br><code>h_pol</code>: <code>True</code> for horizontal polarization, <code>False</code> for vertical polarization </td><td> Specific attenuation (dB/km) </td><td> Forward scattering </td></tr></tbody></table>

<h2>Refractive indices: the <code>refractive</code> module</h2>

This module has various functions to aid with refractive indices, especially in the context of microwave measurements.<br>
<br>
Two methods are available for calculating effective medium approximations (EMA):<br>
<table><thead><th> <b>Method</b> </th><th> <b>Arguments</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> <code>mg_refractive(m, mix)</code> </td><td> <code>m</code>: A tuple giving the refractive indices of the components<br><code>mix</code>: A tuple with the volume fractions of the components (if their sum is not 1, they are normalized with the sum). This should be the same size as <code>m</code>. </td><td> The Maxwell-Garnett effective medium approximation. The arguments <code>m</code> and <code>mix</code> should both have at least two items. If <code>len(m)==2</code>, the first element is taken as the matrix and the second as the inclusion. If <code>len(m)&gt;2</code>, the media are mixed recursively so that the last element is used as the inclusion and the second to last as the matrix, then this mixture is used as the last element on the next iteration, and so on. The effective complex refractive index is then returned. </td></tr>
<tr><td> <code>mg_bruggeman(m, mix)</code> </td><td> <code>m</code>: As above<br><code>mix</code>: As above. </td><td> The Bruggeman effective medium approximation. This works as above, but only supports two components. </td></tr></tbody></table>

There are also definitions of the microwave complex refractive index of water for various wavelengths:<br>
<table><thead><th> <b>Constant</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> <code>m_w_0C</code> </td><td> The complex refractive indices of water in the microwave at 0 °C temperature. This is a dictionary with keys for each of the wavelengths mentioned in the <code>tmatrix_aux</code> description. </td></tr>
<tr><td> <code>m_w_10C</code> </td><td> As above, but for 10 °C temperature. </td></tr>
<tr><td> <code>m_w_20C</code> </td><td> As above, but for 20 °C temperature. </td></tr></tbody></table>

<h1>Examples</h1>

The following examples are adopted from the <code>test.test_tmatrix</code> module. They assume that the following imports have been made:<br>
<pre><code>from pytmatrix.tmatrix import Scatterer<br>
from pytmatrix.psd import PSDIntegrator, GammaPSD<br>
from pytmatrix import orientation, radar, tmatrix_aux, refractive<br>
</code></pre>

A single-orientation case:<br>
<pre><code>&gt;&gt;&gt; scatterer = Scatterer(radius=2.0, wavelength=6.5, m=complex(1.5,0.5), axis_ratio=1.0/0.6)<br>
&gt;&gt;&gt; scatterer.get_SZ()<br>
(array([[  3.89338755e-02 -2.43467777e-01j,<br>
         -1.11474042e-24 -3.75103868e-24j],<br>
       [  1.11461702e-24 +3.75030914e-24j,<br>
         -8.38637654e-02 +3.10409912e-01j]]),<br>
 array([[  8.20899248e-02,  -2.12975199e-02,  -1.94051304e-24,<br>
          2.43057373e-25],<br>
       [ -2.12975199e-02,   8.20899248e-02,   2.00801268e-25,<br>
         -1.07794906e-24],<br>
       [  1.94055633e-24,  -2.01190190e-25,  -7.88399525e-02,<br>
          8.33266362e-03],<br>
       [  2.43215306e-25,  -1.07799010e-24,  -8.33266362e-03,<br>
         -7.88399525e-02]]))<br>
</code></pre>

An orientation-averaged case with a 20° standard deviation Gaussian distribution:<br>
<pre><code>&gt;&gt;&gt; scatterer = Scatterer(radius=2.0, wavelength=6.5, m=complex(1.5,0.5), axis_ratio=1.0/0.6)<br>
&gt;&gt;&gt; scatterer.or_pdf = orientation.gaussian_pdf(std=20.0)<br>
&gt;&gt;&gt; scatterer.orient = orientation.orient_averaged_adaptive<br>
&gt;&gt;&gt; scatterer.get_S()<br>
array([[  6.49005717e-02 -2.42488000e-01j,<br>
         -6.13348029e-16 -4.11094415e-15j],<br>
       [ -1.50045335e-14 -1.63765222e-15j,<br>
         -9.54176591e-02 +2.84758322e-01j]])<br>
</code></pre>

The same, but using fixed-point averaging instead (note the slightly different results):<br>
<pre><code>&gt;&gt;&gt; scatterer = Scatterer(radius=2.0, wavelength=6.5, m=complex(1.5,0.5), axis_ratio=1.0/0.6)<br>
&gt;&gt;&gt; scatterer.or_pdf = orientation.gaussian_pdf(20.0)<br>
&gt;&gt;&gt; scatterer.orient = orientation.orient_averaged_fixed<br>
&gt;&gt;&gt; scatterer.get_S()<br>
array([[  6.49006144e-02 -2.42487917e-01j,<br>
          1.20257280e-11 -5.23021952e-11j],<br>
       [  6.21754954e-12 +2.95662694e-11j,<br>
         -9.54177106e-02 +2.84758152e-01j]])<br>
</code></pre>

A simple example of using particle size distribution integration:<br>
<pre><code>&gt;&gt;&gt; scatterer = Scatterer(wavelength=6.5, m=complex(1.5,0.5), axis_ratio=1.0/0.6)<br>
&gt;&gt;&gt; scatterer.psd_integrator = PSDIntegrator()<br>
&gt;&gt;&gt; scatterer.psd = GammaPSD(D0=1.0, Nw=1e3, mu=4)<br>
&gt;&gt;&gt; scatterer.psd_integrator.D_max = 10.0<br>
&gt;&gt;&gt; scatterer.psd_integrator.init_scatter_table(scatterer)<br>
&gt;&gt;&gt; scatterer.get_Z()<br>
array([[  7.20539942e-02,  -1.54020511e-02,  -9.96222004e-25,<br>
          8.34245116e-26],<br>
       [ -1.54020511e-02,   7.20539942e-02,   1.23279117e-25,<br>
          1.40047875e-25],<br>
       [  9.96224481e-25,  -1.23290932e-25,  -6.89738802e-02,<br>
          1.38873117e-02],<br>
       [  8.34136779e-26,   1.40047656e-25,  -1.38873117e-02,<br>
         -6.89738802e-02]])<br>
</code></pre>

A somewhat more complex example of computing radar scattering properties for a C-band radar, showing more features of the size distribution module as well as radar-specific functions:<br>
<pre><code># For testing variable aspect ratio<br>
# This is from Thurai et al., J. Atmos. Ocean Tech., 24, 2007<br>
def drop_ar(D_eq):<br>
    if D_eq &lt; 0.7:<br>
        return 1.0;<br>
    elif D_eq &lt; 1.5:<br>
        return 1.173 - 0.5165*D_eq + 0.4698*D_eq**2 - 0.1317*D_eq**3 - \<br>
            8.5e-3*D_eq**4<br>
    else:<br>
        return 1.065 - 6.25e-2*D_eq - 3.99e-3*D_eq**2 + 7.66e-4*D_eq**3 - \<br>
            4.095e-5*D_eq**4 <br>
<br>
&gt;&gt;&gt; scatterer = Scatterer(wavelength=tmatrix_aux.wl_C, m=refractive.m_w_10C[tmatrix_aux.wl_C])<br>
&gt;&gt;&gt; scatterer.psd_integrator = PSDIntegrator()<br>
&gt;&gt;&gt; scatterer.psd_integrator.axis_ratio_func = lambda D: 1.0/drop_ar(D)<br>
&gt;&gt;&gt; scatterer.psd_integrator.D_max = 10.0<br>
&gt;&gt;&gt; scatterer.psd_integrator.geometries = (tmatrix_aux.geom_horiz_back, tmatrix_aux.geom_horiz_forw)<br>
&gt;&gt;&gt; scatterer.or_pdf = orientation.gaussian_pdf(20.0)<br>
&gt;&gt;&gt; scatterer.orient = orientation.orient_averaged_fixed<br>
&gt;&gt;&gt; scatterer.psd_integrator.init_scatter_table(scatterer)<br>
&gt;&gt;&gt; scatterer.psd = GammaPSD(D0=2.0, Nw=1e3, mu=4)<br>
&gt;&gt;&gt; radar.refl(scatterer)<br>
6382.9718136764586<br>
&gt;&gt;&gt; radar.refl(scatterer, False)<br>
5066.7211159239041<br>
&gt;&gt;&gt; radar.Zdr(scatterer)<br>
1.2598450941838721<br>
&gt;&gt;&gt; radar.ldr(scatterer)<br>
0.0021935578301794474<br>
&gt;&gt;&gt; radar.rho_hv(scatterer)<br>
0.0021935578301794474<br>
&gt;&gt;&gt; scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)<br>
&gt;&gt;&gt; radar.Kdp(scatterer)<br>
0.1933072695008976<br>
&gt;&gt;&gt; radar.Ai(scatterer)<br>
0.01892301354749585<br>
</code></pre>

<h1>Code documentation</h1>
The public interface of the code is documented with Python docstrings. If you have <a href='http://ipython.org/'>IPython</a>, you can type a question mark after an object to get the docstring. For example:<br>
<pre><code>In [1]: from pytmatrix import radar<br>
<br>
In [2]: radar.Kdp?<br>
Type:       function<br>
Base Class: &lt;type 'function'&gt;<br>
String Form:&lt;function Kdp at 0x40017d0&gt;<br>
Namespace:  Interactive<br>
Definition: radar.Kdp(tm)<br>
Docstring:<br>
Specific differential phase (K_dp) for the current setup.<br>
<br>
Args:<br>
    tm: a TMatrix instance.<br>
<br>
Returns:<br>
   K_dp [deg/km].<br>
<br>
NOTE: This only returns the correct value if the particle diameter and<br>
wavelength are given in [mm]. The tm object should be set to forward<br>
scattering geometry before calling this function.<br>
</code></pre>

<h1>Current limitations and known problems</h1>
<ul><li>The main problem with pytmatrix at the moment is that if an error occurs in the Fortran T-matrix code, this will cause the calculation to exit, <i>also exiting the Python interpreter</i> at the same time. This can happen, e.g., with invalid input values or if the T-matrix code fails to converge. While this does not affect the correctness of the results, it is an inconvenience that one should be aware of.<br>
</li><li>Chebyshev particles are currently not supported as particle shapes (cylinders and spheroids should work normally, but most of the testing has been done with spheroids).</li></ul>

<h1>Citing the code</h1>
No paper publication to cite is currently available. If you'd like to cite this code in a publication, the following reference is suggested:<br>
<br>
Leinonen, J., <i>Python code for T-matrix scattering calculations</i>. Available at <a href='http://code.google.com/p/pytmatrix/'>http://code.google.com/p/pytmatrix/</a>.<br>
<br>
Relevant references should be also made to the Fortran T-matrix code that forms the core of this package. The following should be cited, as appropriate:<br>
<br>
<ol><li>M. I. Mishchenko and L. D. Travis, T-matrix computations of light scattering by large spheroidal particles, <i>Opt. Commun.</i>, vol. 109, 16-21 (1994).<br>
</li><li>M. I. Mishchenko, L. D. Travis, and A. Macke, Scattering of light by polydisperse, randomly oriented, finite circular cylinders, <i>Appl. Opt.</i>, vol. 35, 4927-4940 (1996).<br>
</li><li>D. J. Wielaard, M. I. Mishchenko, A. Macke, and B. E. Carlson,       Improved T-matrix computations for large, nonabsorbing and weakly absorbing nonspherical particles and comparison with geometrical optics approximation, <i>Appl. Opt.</i>, vol. 36, 4305-4313 (1997).</li></ol>

<h1>Original code</h1>
The original Fortran 77 code is available at: <a href='http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html'>http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html</a>.<br>
<br>
<h1>References</h1>
The following books were used extensively as source material to build this code:<br>
<br>
M. I. Mishchenko, J. W. Hovenier, and L. D. Travis (ed.), <i>Light Scattering by Nonspherical Particles</i>, Academic Press (2000).<br>
<br>
Bringi, V. N., and V. Chandrasekar, <i>Polarimetric Doppler weather radar: principles and applications</i>, Cambridge University Press (2001).<br>
<br>
A general review of the T-matrix approach can be found in:<br>
<br>
M. I. Mishchenko, L. D. Travis, and D. W. Mackowski, T-matrix computations of light scattering by nonspherical particles:  a review, <i>J. Quant. Spectrosc. Radiat. Transfer</i>, vol. 55, 535-575 (1996).<br>
<br>
Additional useful information is contained in the paper:<br>
<br>
M. I. Mishchenko and L. D. Travis, Capabilities and limitations of a current FORTRAN implementation of the T-matrix method for randomly oriented, rotationally symmetric scatterers, <i>J. Quant. Spectrosc. Radiat. Transfer</i>, vol. 60, 309-324 (1998).<br>
<br>
The definitions and notation used can also be found in:<br>
<br>
M. I. Mishchenko, Calculation of the amplitude matrix for a nonspherical particle in a fixed orientation, <i>Appl. Opt.</i>, vol. 39, 1026-1031 (2000).<br>
<br>
<h1>Version history</h1>

0.2:<br>
<ul><li>Added calculation of generic scattering parameters (scattering cross section, extinction cross section, single-scattering albedo, asymmetry parameter).<br>
</li><li>Changed how PSD integration works (see documentation above).<br>
</li><li>Replaced quadrature-point code with a pure Python version.</li></ul>

0.1: First release<br>
<br>
<h1>See also</h1>
Python code for calculating Mie scattering from single- and dual-layered spheres: <a href='http://code.google.com/p/pymiecoated/'>pymiecoated</a>.
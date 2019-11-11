#### Guide to installing MMSEQ dependencies

- [Boost C++ Libraries](http://boost.org/)
  - Mac OS X+[MacPorts](http://www.macports.org/): `port install boost`
  - Mac OS X+[Homebrew](http://brew.sh/): `brew install boost; ln -s /usr/local/lib/libboost_regex-mt.a /usr/local/lib/libboost_regex.a; ln -s /usr/local/lib/libboost_iostreams-mt.a /usr/local/lib/libboost_iostreams.a`
  - Ubuntu/Debian/etc: `apt-get install libboost-all-dev`
  - Fedora/CentOS/Red Hat/etc: `yum install boost` (or, on older systems, `yum install boost141`)
- [GNU Scientific Library](http://www.gnu.org/software/gsl)
  - Mac OS X+[MacPorts](http://www.macports.org/): `port install gsl`
  - Mac OS X+[Homebrew](http://brew.sh/): `brew install gsl`
  - Ubuntu/Debian/etc: `apt-get install libgsl0-dev`
  - Fedora/CentOS/Red Hat/etc: `yum install gsl`
- [Armadillo C++ linear algebra library](http://arma.sf.net)
  - Mac OS X+[Homebrew](http://brew.sh/): `brew install homebrew/science/openblas; brew install homebrew/dupes/lapack; brew install armadillo`
- [HTSlib library](http://htslib.org/)
  - Mac OS X+[Homebrew](http://brew.sh/): `brew install htslib`

Note that due to a lack of OpenMP support in Apple's clang compiler (as of El Capitan), the Mac binaries will be single-threaded.

Make sure you have a [Ruby interpreter](http://ruby-lang.org/) to run the Ruby scripts:

- Ubuntu/Debian/etc: `apt-get install ruby-full`
- Fedora/CentOS/Red Hat/etc: `yum install ruby`

To correct the transcript lengths for non-uniformity effects using the Poisson regression method of [Li et al.](http://genomebiology.com/2010/11/5/R50), install the R package [mseq](http://cran.r-project.org/web/packages/mseq/index.html).

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
- [SAMtools library](http://samtools.sourceforge.net/)
  - Mac OS X+[Homebrew](http://brew.sh/): `brew install samtools; export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include/bam`

Due to a bug in the default version of GCC that comes with XCode for Mac OS X, you need to use [a more recent version of GCC](http://hpc.sourceforge.net/). If libraries aren't found by the compiler, add the path(s) containing the header files and libraries respectively to the `LIBRARY_PATH` and `CPLUS_INCLUDE_PATH` environment variables. E.g. to add the path to SAMtools (which should contain `sam.h` and `libbam.a`, among other files):

    export CPLUS_INCLUDE_PATH=/home/bob/seq/samtools-0.1.17:$CPLUS_INCLUDE_PATH
    export LIBRARY_PATH=/home/bob/seq/samtools-0.1.17:$LIBRARY_PATH

Also make sure you have a [Ruby interpreter](http://ruby-lang.org/) to run the Ruby scripts:

- Ubuntu/Debian/etc: `apt-get install ruby-full`
- Fedora/CentOS/Red Hat/etc: `yum install ruby`

To correct the transcript lengths for non-uniformity effects using the Poisson regression method of [Li et al.](http://genomebiology.com/2010/11/5/R50), install the R package [mseq](http://cran.r-project.org/web/packages/mseq/index.html).

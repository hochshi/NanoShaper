FROM rockylinux:9
RUN echo "Installing Nanoshaper 1.5"

RUN dnf install -y boost-devel
RUN dnf install -y boost
RUN dnf install -y gmp
RUN dnf install -y gmp-devel
RUN dnf install -y mpfr
RUN dnf install -y mpfr-devel
RUN dnf install -y gcc
RUN dnf install -y gcc-c++
RUN dnf install -y make
RUN dnf install -y tar
RUN dnf install -y git
RUN dnf install -y wget
RUN dnf install -y unzip
RUN dnf install -y openssl-devel
RUN dnf install -y cmake

###TBB###
WORKDIR /opt/
RUN wget https://github.com/wjakob/tbb/archive/refs/heads/master.zip
RUN unzip master.zip
RUN mv tbb-master tbb
WORKDIR /opt/tbb
RUN rm -f CMakeCache.txt
RUN cmake .
RUN make -j8
RUN make install
WORKDIR /opt/
RUN rm master.zip

###CGAL-5.6.2###
ENV CGALVER="5.6.2"
ENV CGALFILE="v5.6.2.tar.gz"
ENV CGALDIR="cgal-5.6.2"
WORKDIR /opt/
RUN wget https://github.com/CGAL/cgal/archive/v${CGALVER}.tar.gz
RUN tar xvf ${CGALFILE}
RUN rm ${CGALFILE}
WORKDIR /opt/${CGALDIR}
RUN cmake . -DBUILD_SHARED_LIBS=OFF -DWITH_examples=false -DWITH_CGAL_Qt5=false -DWITH_CGAL_Qt4=false -DWITH_CGAL_Qt3=false -DWITH_CGAL_ImageIO=false 
RUN make clean
RUN make

###NanoShaper###
RUN mkdir /usr/local/nanoshaper/
WORKDIR /usr/local/nanoshaper/
RUN git clone https://gitlab.iit.it/SDecherchi/nanoshaper.git /usr/local/nanoshaper/
RUN git checkout NanoShaper_PatchBased
RUN rm -fr example
RUN cp CMakeLists_standalone.txt CMakeLists.txt 
WORKDIR /usr/local/nanoshaper/build
RUN cmake .. -DCGAL_DIR=/opt/${CGALDIR} -DCMAKE_BUILD_TYPE="Release"
RUN make clean 
RUN make -j8

# a bit compressed image with only the NS executable
FROM rockylinux:9
RUN echo "Release Image - Nanoshaper 1.5"
RUN dnf update -y
RUN dnf install -y boost-devel
RUN dnf install -y boost
RUN dnf install -y gmp
RUN dnf install -y gmp-devel
RUN dnf install -y mpfr
RUN dnf install -y mpfr-devel
RUN mkdir /usr/local/nanoshaper
RUN mkdir /usr/local/nanoshaper/build
WORKDIR /usr/local/nanoshaper/build
COPY --from=0 /usr/local/nanoshaper/build/NanoShaper .
COPY --from=0 /usr/local/lib/libtbb.so /usr/local/lib/
COPY --from=0 /usr/local/lib/libtbbmalloc.so /usr/local/lib/
COPY --from=0 /usr/local/lib/libtbbmalloc_proxy.so /usr/local/lib/


ENV CONFFILE="conf.prm"

VOLUME ["/App"]
WORKDIR /App
ENTRYPOINT ["/usr/local/nanoshaper/build/NanoShaper"]
CMD [${CONFFILE}]

FROM thomasbigot/bpp:v2.2.0
MAINTAINER Thomas Bigot <tomed.bigot@gmail.com>
RUN git clone git://pbil-dev.univ-lyon1.fr/phyldog
RUN mkdir phyldog/build
RUN apt-get install -y cmake
RUN apt-get install -y g++-4.9
RUN /bin/bash -c 'wget http://www.libpll.org/Downloads/libpll-1.0.11.tar.gz ;\
tar xvzf libpll-1.0.11.tar.gz;\
cd libpll-1.0.11;\
export CXX=/usr/bin/g++-4.9;\
./configure;\
make install;\'
RUN apt-get -y install libboost-all-dev
RUN cd phyldog/build && cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9
RUN cd phyldog/build && make install 

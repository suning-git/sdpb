This guide is for building SDPD.

* [Installation](#installation)

# Installation

1. sdpd is based on sdpb 2.4.0. Try to built sdpd 2.4.0 first and all its required libraries.

2. On the cluster of Perimeter Institute, I use following commands to build sdpd

        module load openmpi/gcc-9/64/4.0.5
        module load cmake/3.18.1
        module load boost/1.67.0
        module load python/3.7
        module load trilinos/13.0.0
        module load openblas
        
        git clone https://github.com/suning1985/sdpb.git
        cd sdpb
        git checkout sdpd_experimental
        
        ./waf configure --prefix /home/nsu2/packages/install/sdpd --elemental-dir=/cm/shared/apps/elemental/sdpb-fork/ --boost-dir $BOOST_ROOT --rapidjson-incdir=/cm/shared/apps/rapidjson/1.1.0/include
        ./waf install

    The commands are essentially the same from sdpb 2.4.0, except downloading. A binary named sdpd will appear in the folder specified by --prefix.

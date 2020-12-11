# to set env for esmf/modules
# source this file: esmf-impi-env.sh
#

# Path to libraries and includes and bins
module load intel/18.0.5.274 
module load szip/2.1  
module load hdf5/1.10.4 
module load impi/2018.4.274
module load netcdf/4.6.1


# Environment for ESMF v8.0.0 beta snapshot 40
module use /home/emc.nemspara/SOFT-hera/modulefiles
# module load intel/15.1.133 impi/5.1.1.109 netcdf/4.3.0 
# module load impi/5.1.1.109
#module load yaml-cpp 
module load esmf/8.0.0bs48g

# latest version on Here
# module load yaml-cpp
# module --ignore-cache load "esmf/8.0.0bs42"
# module load esmf/7.1.0r

# Environment for ESMF 
# export ESMF_DIR=/scratch2/COASTAL/coastal/save/Beheen.M.Trimble/esmf_v8.0_beta/INSTALL/esmf8.0.beta_impi18.0.4
# export ESMF_BOPT='g'
# export ESMF_COMM=intelmpi      # mpich, mpich2,lam, openmpi or intelmpi
# export ESMF_COMPILER='intel'
# export ESMF_ABI='64'
# export ESMF_OS='Linux'         # uname -s

# export ESMFMKFILE=$ESMF_DIR/lib/esmf.mk

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch2/COASTAL/coastal/save/Beheen.M.Trimble/yaml-cpp-0.6.2/lib:$ESMF_DIR/lib


## Contents

* [Usage](#installation-and-usage)

# SDPD

SDPD computes `dx,dy,dX,dY` and derivative of the primal objective with respect to the change of bootstrap condition.

## Usage

The usage of sdpd is

    sdpd --procsPerNode=[N] --[MODE] -s [SDP1] -d [SDP2] -i [INCK] -o [OUTCK] --precision=[PRECISION]

`[N]` is the number of cpu cores. `[MODE]` specify which mode of sdpd, see below. `[SDP1]` and `[SDP2]` are two sdp folders. `[INCK]` is a checkpoint folder contains x,y,X,Y in text format (from sdpb `--writeSolution="x,y,X,Y"` option). `[OUTCK]` is a folder for save `dx,dy,dX,dY`. The `-o` is optional. If not specified, `dx,dy,dX,dY` will not be saved. `[PRECISION]` should be the same precision that you used in sdpb to get `[INCK]`.

Before above command, certain MPI command should be specified. The full command might look like `mpirun -n [N] sdpd...` . The `mpirun -n [N]` part could be slightly different in different clusters. It should be the command you use for sdpb 2.4.0.

SDPD can work in 3 different modes. In all modes, the derivative of the primal objective is returned by printing out `[SDPDReturnBegin]...[SDPDReturnEnd]`.

### Mode 1: `--SDP1_db`

In this mode, `[SDP1]` is a sdp folder contains same bootstrap condition that used to generate `[INCK]`, except the objective `b` is replace by a derivative `db`. `-d` should not be specified. SDPD will compute the `dx,dy,dX,dY` with respect to change of `db`.

Example :

    mpirun -n 80 sdpd --procsPerNode=80 --SDP1_db -s ./Proj_SDPDExample/SDPBFiles/sdp1db.sdp -i ./Proj_SDPDExample/SDPBFiles/sdp1theta.out --precision 768


### Mode 2: `--SDP2_B_b_c`

In this mode, `[SDP1]` is a sdp folder contains same bootstrap condition that used to generate `[INCK]`. `[SDP2]` is a sdp folder contains slightly deformed bootstrap condition. SDPD will subtract SDP1-SDP2 to obtain `dB,db,dc` and compute the `dx,dy,dX,dY` with respect to those changes.

Example :

    mpirun -n 80 sdpd --procsPerNode=80 --SDP2_B_b _c -s ./Proj_SDPDExample/SDPBFiles/sdp1theta.sdp -d ./Proj_SDPDExample/SDPBFiles/sdp2theta.sdp -i ./Proj_SDPDExample/SDPBFiles/sdp1theta.out -o ./Proj_SDPDExample/SDPBFiles/sdp1theta.sdpd.out --precision 768


### Mode 3: `--SDP2_dB_db_dc`

In this mode, `[SDP1]` is a sdp folder contains same bootstrap condition that used to generate `[INCK]`. `[SDP2]` is a sdp folder contains the change of bootstrap condition `dB,db,dc`. SDPD computes the `dx,dy,dX,dY` with respect to those changes.

Example :

    mpirun -n 80 sdpd --procsPerNode=80 --SDP2_dB_db _dc -s ./Proj_SDPDExample/SDPBFiles/sdp1theta.sdp -d ./Proj_SDPDExample/SDPBFiles/dsdptheta.sdp -i ./Proj_SDPDExample/SDPBFiles/sdp1theta.out --precision 768


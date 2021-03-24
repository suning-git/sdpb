The test sdp problem is parametrized by two parameter s,e. Side note : The test sdp problem is from 2D Ising single correlator bootstrap, where s,e are scale dimension of \simga and \epsilon.

### sdpb run

In this folder, sdp1.xml, sdp2.xml, sdp3.xml are XMLs for sdp1=sdp(s0,e0), sdp2=sdp(s0-10^-20,e0), sdp(s0,e0-10^-20) respectively (side note : s0=0.125, e0=1). Between sdp1 and sdp2, the SDP variable B,b,c changes. Between sdp1 and sdp3, only B changes.

Convert the 3 xml to sdp using pvm2sdp, then run they with following sdpb parameters :

    --dualityGapThreshold=1e-30 --primalErrorThreshold=1e-20 --dualErrorThreshold=1e-20 --precision=768 --initialMatrixScalePrimal=1e+20 --initialMatrixScaleDual=1e+20 --maxComplementarity=1e+100 --writeSolution="x,y,X,Y"

We found the primal objectives for the XMLs are:

sdp1.xml :
-6.5908836257341204303178949252821185459109323930172307657699609106454761117100204394775217644928635671684301496655151126378116169095945301434599395810510934599359016901363864212674348204098751672147817938395796996079410919715026084439*10^-6

sdp2.xml :
-6.5908836257341204293049488830315420296612669551331319703473620248886579323519405230404295815904065651511428261972351097276728196722891230000000000000000000000000000000000000000000000000000000000000000*10^-6

sdp3.xml :
-6.5908836257341204293049488830315420296612669551331319703473620248886579323519405230404295815904065651511428261972351097276728196722891230000000000000000000000000000000000000000000000000000000000000000*10^-6


### sdp_derivative mode 2 `--SDP2_B_b_c` run

Using following command to run sdp_derivative in `--SDP2_B_b_c` mode:

    mpirun -n 80 sdpd --procsPerNode=80 --SDP2_B_b_c -s ./sdp1.sdp -d ./sdp2.sdp -i ./sdp1.out --precision 768
    mpirun -n 80 sdpd --procsPerNode=80 --SDP2_B_b_c -s ./sdp1.sdp -d ./sdp3.sdp -i ./sdp1.out --precision 768

Extract the number from [SDPDReturnBegin]...[SDPDReturnEnd], we found the sdp_derivative mode 2 for the derivative of the primal objectives are
In s direction (using sdp2) : -1.0129460422505765115066347412619159636036371724689919343195150436621308212350723362474714447574304974753225690093385135517725516934321577213194773779860591443441891190236154268531745702432151896346338214461986734551851787364261472866*10^-24
In e direction (using sdp3) : -2.536074828321516824601161829520163417793265719403734172392866708779501471320438232767914724758626907757228070123929252423369646933105633836635557541482874801451464849414870459281213083528136394867702684769935585732070837092572039352*10^-25

Rescaling those numbers by 10^-20, we found they match with numeric derivative from previous section

### sdp_derivative mode 1 `--SDP1_db` run

The sdp1db.xml contains exact derivative db for d(sdp(s,e))/d(e). Using pvm2sdp to merge this xml with sdp1.xml, i.e. replace the objective in sdp1.xml. Run sdp_derivative with following parameter

    mpirun -n 80 sdpd --procsPerNode=80 --SDP1_db -s ./sdp1db.sdp -i ./sdp1.out --precision 768

Extract the number from [SDPDReturnBegin]...[SDPDReturnEnd], we found the sdp_derivative mode 1 result for the derivative of the primal objectives are

-0.00002536074828321516824651322166944606741931475322377681999859872108622465849301132384020794355242040875994485419243395575754110827749767198791990416582632479397357303647737262991291937698967288987217823596262178443435430549447795700128

It matches with the numerical derivative and mode 2 derivative.

### sdp_derivative mode 3 `--SDP2_dB_db_dc` run

The dsdp.xml contains exact derivative dB,db,dc for d(sdp(s,e))/d(s). Using pvm2sdp to convert this xml to dsdp.sdp and run sdp_derivative with following parameter

    mpirun -n 80 sdpd --procsPerNode=80 --SDP2_dB_db_dc -s ./sdp1.sdp -d ./dsdp.sdp -i ./sdp1.out --precision 768
    
Extract the number from [SDPDReturnBegin]...[SDPDReturnEnd], we found the sdp_derivative mode 3 result for the derivative of the primal objectives are

-0.00010129460422505765163316139840277616441596306857897382196559332253399035234445073699960597620353893890898447973557947725298971751384869807193530715617873265205937602127505955910091901283177380365472819614827427761649678488889820943354

It matches with the numerical derivative and mode 2 derivative.

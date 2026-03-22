#ifndef RTW_HEADER_sistema_completo_h_
#define RTW_HEADER_sistema_completo_h_
#ifndef sistema_completo_COMMON_INCLUDES_
#define sistema_completo_COMMON_INCLUDES_
#include <stdlib.h>
#include "sl_AsyncioQueue/AsyncioQueueCAPI.h"
#include "rtwtypes.h"
#include "sigstream_rtw.h"
#include "simtarget/slSimTgtSigstreamRTW.h"
#include "simtarget/slSimTgtSlioCoreRTW.h"
#include "simtarget/slSimTgtSlioClientsRTW.h"
#include "simtarget/slSimTgtSlioSdiRTW.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "raccel.h"
#include "slsv_diagnostic_codegen_c_api.h"
#include "rt_logging_simtarget.h"
#include "dt_info.h"
#include "ext_work.h"
#endif
#include "sistema_completo_types.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include <stddef.h>
#include "rtGetInf.h"
#include "rtw_modelmap_simtarget.h"
#include "rt_defines.h"
#include <string.h>
#define MODEL_NAME sistema_completo
#define NSAMPLE_TIMES (3) 
#define NINPUTS (0)       
#define NOUTPUTS (0)     
#define NBLOCKIO (129) 
#define NUM_ZC_EVENTS (0) 
#ifndef NCSTATES
#define NCSTATES (12)   
#elif NCSTATES != 12
#error Invalid specification of NCSTATES defined in compiler command
#endif
#ifndef rtmGetDataMapInfo
#define rtmGetDataMapInfo(rtm) (*rt_dataMapInfoPtr)
#endif
#ifndef rtmSetDataMapInfo
#define rtmSetDataMapInfo(rtm, val) (rt_dataMapInfoPtr = &val)
#endif
#ifndef IN_RACCEL_MAIN
#endif
typedef struct { real_T k2sn1moatl ; real_T m4zsvrd15a ; real_T jpgl5uqnen ;
real_T hpaqwswri4 ; real_T ojegsw1r4w ; real_T hpqeqr0o0i ; real_T efmorfp5hk
; real_T homtptr1wn ; real_T dbz2he40vm ; real_T mbfoowmnkk ; real_T
dlox32o0ua ; real_T oh5qhuuk4l ; real_T dx12nyer4s ; real_T c5fdscj3uw ;
real_T lnhoszrvh2 ; real_T c11k2clsrb ; real_T mqiqfyapeu ; real_T kylrj4an0z
; real_T dshebsf5x3 ; real_T ly2xklr4fz ; real_T h1lgwz5ovu ; real_T
gifca5zh5d ; real_T gysjfnsk1k ; real_T ajowt5z1q5 ; real_T fwvyoitwh3 ;
real_T netv40pb4f ; real_T jo0ai5lhtc ; real_T k532jfpdyn ; real_T nxola2u2g3
; real_T kmls4v4rl0 ; real_T af03zirre2 ; real_T dhufv0kct4 ; real_T
fd152jijlt ; real_T klsfhtexlj ; real_T p3yvsgsymc ; real_T elrtlardao ;
real_T pxphckphjm ; real_T k25kf23xp1 ; real_T g5whoifsex ; real_T fb0zj32vbu
; real_T kga0ttwxih ; real_T kf3sjuqaae ; real_T fwd4ryyrfq [ 3 ] ; real_T
iamwwft4eh ; real_T edszxpjalm ; real_T azzpx0o04y ; real_T nd4l40ce0g ;
real_T pr2t5fqsxs ; real_T nly0v00iut ; real_T cri0rzhvjd ; real_T o1eaam1gvh
; real_T i1cztljifv ; real_T hst2mc5opt ; real_T kpwiznui1d ; real_T
ollu0ludgl ; real_T bpf3qazr44 ; real_T nc1t1u3huh ; real_T dnzl4cfy5s ;
real_T ougga4l4aa ; real_T pxdilwla4w ; real_T gl1fdweeu4 ; real_T fhjbzzgwbv
; real_T kl3fuewomp ; real_T h1yubcelv5 ; real_T k1afv5mxmb ; real_T
lggf4lsbsd ; real_T evhxefbhrl ; real_T gt3nq2cscb ; real_T aqulfpaq14 ;
real_T bd4jar1csl ; real_T dvu22fpphc ; real_T fnvmvvqodz [ 3 ] ; real_T
dz4hiszfdk ; real_T lswxnkpgli ; real_T oqeh2ikwco ; real_T jxvum3zvoh ;
real_T ou4cqm5an5 ; real_T onoapnnemw ; real_T iktm1emn42 ; real_T kwquioc0a2
; real_T caar3rhznz ; real_T gbixzpfrni ; real_T jwtzlohzhs ; real_T
cbv3kpprwx ; real_T cdjjap5tb4 ; real_T k5iqwa5dub ; real_T lbeji2oewc ;
real_T crxm0t2udf ; real_T avres0c4zw ; real_T buftnulajc ; real_T pwtikpf0j0
; real_T marzp3jauv ; real_T dymq05epmd ; real_T d3woaso1vk ; real_T
lata252i34 ; real_T p5w0nljbtt ; real_T jkeaxbjqzp ; real_T lpe0g0uuxy ;
real_T kwzodnsn51 ; real_T kjg2pniu00 ; real_T c4ogmaamk3 ; real_T jivucflgxz
; real_T b0djz0xuyx ; real_T moi1nwsely ; real_T apbcudoiru ; real_T
bptosmp44x ; real_T fxltoioz2o ; real_T k3waxxng0u ; real_T ek2tbc2vjr ;
real_T j1emwd2qae ; real_T lsyqe5i55w ; real_T c4n4jb2dgo ; real_T ny2zw5nsgp
; real_T fzxuttsywh ; real_T lb5025wpa3 ; real_T ib34wfcjm3 ; real_T
nf03pyptgh ; real_T iss1dct3mr ; real_T hpveetbx1j ; real_T eufxrbxe4z ;
real_T j0wuzgingg ; real_T aduksi4izn ; real_T pia0trm2jv ; } B ; typedef
struct { real_T parg3eaod2 ; real_T na0wt4oepn ; real_T dejwxn1xy1 ; real_T
dcb22robm1 ; struct { void * TimePtr ; void * DataPtr ; void * RSimInfoPtr ;
} jrzquc22ff ; struct { void * TimePtr ; void * DataPtr ; void * RSimInfoPtr
; } ml4gzok05f ; struct { void * TimePtr ; void * DataPtr ; void *
RSimInfoPtr ; } diriss0ycm ; struct { void * LoggedData ; } dsna2h3s4j ;
struct { void * LoggedData ; } pv5y2snvjf ; struct { void * LoggedData ; }
imrctj1fis ; struct { void * LoggedData ; } ak1r52kzjl ; struct { void *
LoggedData ; } o0dx25v4bj ; struct { void * AQHandles ; } opv5nm4cbf ; struct
{ void * AQHandles ; } cqxqnkszbm ; struct { void * AQHandles ; } dt1kwx4eay
; struct { void * AQHandles ; } dnhkap2j1u ; struct { void * LoggedData ; }
put0eavhrg ; struct { void * LoggedData ; } n02okhcqno ; struct { void *
LoggedData [ 3 ] ; } a0ojjrdjng ; struct { void * LoggedData [ 3 ] ; }
ebooqd5zoc ; struct { void * LoggedData ; } i2e1lohrez ; struct { void *
LoggedData ; } mvlzl14qts ; struct { void * TimePtr ; void * DataPtr ; void *
RSimInfoPtr ; } ekcesp1k2p ; struct { void * TimePtr ; void * DataPtr ; void
* RSimInfoPtr ; } hrwve4sncv ; struct { void * AQHandles ; } pqdji2nimh ;
int32_T asp1sby2rd ; int32_T ekogk2wdou ; int32_T axuxe1t12v ; int32_T
fpyg5ifl2j ; struct { int_T PrevIndex ; } jptvtcgs5c ; struct { int_T
PrevIndex ; } l2ca51dh03 ; struct { int_T PrevIndex ; } h33gajte4x ; struct {
int_T PrevIndex ; } gjvhwykizl ; struct { int_T PrevIndex ; } ojtfq2nxq4 ;
uint8_T iqsjj14to5 ; uint8_T derxydwg12 ; uint8_T ef11rieqzy ; uint8_T
aa1tc5fl4n ; boolean_T hpxndlhxow ; boolean_T d3jkgjc1oh ; boolean_T
gh0bm4uzee ; boolean_T lljtkeuzy4 ; } DW ; typedef struct { real_T o5hpeumlsb
; real_T jj3u5nl5xw ; real_T gjk4wr4fl2 ; real_T mmd3jbgtkz ; real_T
lwhlamef2l ; real_T logvqwtjje ; real_T mwqg0fpixr ; real_T fj4qnnj32q ;
real_T kup1qcjf04 ; real_T asvy4gz4nf ; real_T ncrlhqqtwv ; real_T b2dj4ozsrg
; } X ; typedef struct { real_T o5hpeumlsb ; real_T jj3u5nl5xw ; real_T
gjk4wr4fl2 ; real_T mmd3jbgtkz ; real_T lwhlamef2l ; real_T logvqwtjje ;
real_T mwqg0fpixr ; real_T fj4qnnj32q ; real_T kup1qcjf04 ; real_T asvy4gz4nf
; real_T ncrlhqqtwv ; real_T b2dj4ozsrg ; } XDot ; typedef struct { boolean_T
o5hpeumlsb ; boolean_T jj3u5nl5xw ; boolean_T gjk4wr4fl2 ; boolean_T
mmd3jbgtkz ; boolean_T lwhlamef2l ; boolean_T logvqwtjje ; boolean_T
mwqg0fpixr ; boolean_T fj4qnnj32q ; boolean_T kup1qcjf04 ; boolean_T
asvy4gz4nf ; boolean_T ncrlhqqtwv ; boolean_T b2dj4ozsrg ; } XDis ; typedef
struct { real_T o5hpeumlsb ; real_T jj3u5nl5xw ; real_T gjk4wr4fl2 ; real_T
mmd3jbgtkz ; real_T lwhlamef2l ; real_T logvqwtjje ; real_T mwqg0fpixr ;
real_T fj4qnnj32q ; real_T kup1qcjf04 ; real_T asvy4gz4nf ; real_T ncrlhqqtwv
; real_T b2dj4ozsrg ; } CStateAbsTol ; typedef struct { real_T o5hpeumlsb ;
real_T jj3u5nl5xw ; real_T gjk4wr4fl2 ; real_T mmd3jbgtkz ; real_T lwhlamef2l
; real_T logvqwtjje ; real_T mwqg0fpixr ; real_T fj4qnnj32q ; real_T
kup1qcjf04 ; real_T asvy4gz4nf ; real_T ncrlhqqtwv ; real_T b2dj4ozsrg ; }
CXPtMin ; typedef struct { real_T o5hpeumlsb ; real_T jj3u5nl5xw ; real_T
gjk4wr4fl2 ; real_T mmd3jbgtkz ; real_T lwhlamef2l ; real_T logvqwtjje ;
real_T mwqg0fpixr ; real_T fj4qnnj32q ; real_T kup1qcjf04 ; real_T asvy4gz4nf
; real_T ncrlhqqtwv ; real_T b2dj4ozsrg ; } CXPtMax ; typedef struct {
rtwCAPI_ModelMappingInfo mmi ; } DataMapInfo ; struct P_ { real_T C_ts ;
real_T J_eq ; real_T K_sa ; real_T K_sia ; real_T L_d ; real_T L_ls ; real_T
L_q ; real_T Pp ; real_T R_s_ref ; real_T R_ts ; real_T T_s_ref ; real_T
alpha_Cu ; real_T b_a ; real_T b_eq ; real_T k_i ; real_T k_l ; real_T
k_omega2 ; real_T k_theta2 ; real_T lambda_m ; real_T r ; real_T
fromWS_Signal1_Time0 [ 2 ] ; real_T fromWS_Signal1_Data0 [ 2 ] ; real_T
Integrator_IC ; real_T Integrator1_IC ; real_T Integrator_IC_i2wxdirj0b ;
real_T Integrator1_IC_oww15y3mqv ; real_T Integrator1_IC_gd0uumlv2e ; real_T
Gain4_Gain ; real_T FromWorkspace_Time0 [ 141 ] ; real_T FromWorkspace_Data0
[ 141 ] ; real_T FromWorkspace_Time0_o2ayzzhu4m [ 2 ] ; real_T
FromWorkspace_Data0_knbhfzr1k4 [ 2 ] ; real_T Integrator_IC_gpeemvg0a0 ;
real_T Integrator1_IC_epj0vehlxx ; real_T Integrator2_IC ; real_T Gain_Gain ;
real_T Gain1_Gain ; real_T Gain2_Gain ; real_T Integrator2_IC_cltf442w1a ;
real_T Gain3_Gain ; real_T Integrator2_IC_lurynlyjlu ; real_T
Integrator_IC_n5atnjvnfo ; real_T Gain_Gain_bpfw2sk42p ; real_T
Gain1_Gain_o051fegvds ; real_T Gain2_Gain_fpelshvuq5 ; real_T
FromWorkspace_Time0_bh3fvyhkqh [ 2 ] ; real_T FromWorkspace_Data0_ahxyd5aapb
[ 2 ] ; real_T Gain1_Gain_bfommv4lak ; real_T Gain2_Gain_gt5xi45qct ; real_T
FromWorkspace_Time0_fivuiqtoah [ 2 ] ; real_T FromWorkspace_Data0_kevgbgk3mx
[ 2 ] ; real_T Integrator2_IC_dzcpbao03e ; real_T Constant_Value ; real_T
Constant2_Value ; uint8_T ManualSwitch1_CurrentSetting ; uint8_T
ManualSwitch_CurrentSetting ; } ; extern const char_T *
RT_MEMORY_ALLOCATION_ERROR ; extern B rtB ; extern X rtX ; extern DW rtDW ;
extern P rtP ; extern mxArray * mr_sistema_completo_GetDWork ( ) ; extern
void mr_sistema_completo_SetDWork ( const mxArray * ssDW ) ; extern mxArray *
mr_sistema_completo_GetSimStateDisallowedBlocks ( ) ; extern const
rtwCAPI_ModelMappingStaticInfo * sistema_completo_GetCAPIStaticMap ( void ) ;
extern SimStruct * const rtS ; extern const int_T gblNumToFiles ; extern
const int_T gblNumFrFiles ; extern const int_T gblNumFrWksBlocks ; extern
rtInportTUtable * gblInportTUtables ; extern const char * gblInportFileName ;
extern const int_T gblNumRootInportBlks ; extern const int_T
gblNumModelInputs ; extern const int_T gblInportDataTypeIdx [ ] ; extern
const int_T gblInportDims [ ] ; extern const int_T gblInportComplex [ ] ;
extern const int_T gblInportInterpoFlag [ ] ; extern const int_T
gblInportContinuous [ ] ; extern const int_T gblParameterTuningTid ; extern
DataMapInfo * rt_dataMapInfoPtr ; extern rtwCAPI_ModelMappingInfo *
rt_modelMapInfoPtr ; void MdlOutputs ( int_T tid ) ; void
MdlOutputsParameterSampleTime ( int_T tid ) ; void MdlUpdate ( int_T tid ) ;
void MdlTerminate ( void ) ; void MdlInitializeSizes ( void ) ; void
MdlInitializeSampleTimes ( void ) ; SimStruct * raccel_register_model (
ssExecutionInfo * executionInfo ) ;
#endif

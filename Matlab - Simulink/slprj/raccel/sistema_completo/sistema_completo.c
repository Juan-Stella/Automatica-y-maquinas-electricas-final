#include "sistema_completo.h"
#include "rtwtypes.h"
#include "mwmathutil.h"
#include "sistema_completo_private.h"
#include "rt_logging_mmi.h"
#include "sistema_completo_capi.h"
#include "sistema_completo_dt.h"
extern void * CreateDiagnosticAsVoidPtr_wrapper ( const char * id , int nargs
, ... ) ; RTWExtModeInfo * gblRTWExtModeInfo = NULL ; void
raccelForceExtModeShutdown ( boolean_T extModeStartPktReceived ) { if ( !
extModeStartPktReceived ) { boolean_T stopRequested = false ;
rtExtModeWaitForStartPkt ( gblRTWExtModeInfo , 2 , & stopRequested ) ; }
rtExtModeShutdown ( 2 ) ; }
#include "slsv_diagnostic_codegen_c_api.h"
#include "slsa_sim_engine.h"
const int_T gblNumToFiles = 0 ; const int_T gblNumFrFiles = 0 ; const int_T
gblNumFrWksBlocks = 5 ;
#ifdef RSIM_WITH_SOLVER_MULTITASKING
boolean_T gbl_raccel_isMultitasking = 1 ;
#else
boolean_T gbl_raccel_isMultitasking = 0 ;
#endif
boolean_T gbl_raccel_tid01eq = 0 ; int_T gbl_raccel_NumST = 3 ; const char_T
* gbl_raccel_Version = "9.9 (R2023a) 19-Nov-2022" ; void
raccel_setup_MMIStateLog ( SimStruct * S ) {
#ifdef UseMMIDataLogging
rt_FillStateSigInfoFromMMI ( ssGetRTWLogInfo ( S ) , & ssGetErrorStatus ( S )
) ;
#else
UNUSED_PARAMETER ( S ) ;
#endif
} static DataMapInfo rt_dataMapInfo ; DataMapInfo * rt_dataMapInfoPtr = &
rt_dataMapInfo ; rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr = & (
rt_dataMapInfo . mmi ) ; const int_T gblNumRootInportBlks = 0 ; const int_T
gblNumModelInputs = 0 ; extern const char * gblInportFileName ; extern
rtInportTUtable * gblInportTUtables ; const int_T gblInportDataTypeIdx [ ] =
{ - 1 } ; const int_T gblInportDims [ ] = { - 1 } ; const int_T
gblInportComplex [ ] = { - 1 } ; const int_T gblInportInterpoFlag [ ] = { - 1
} ; const int_T gblInportContinuous [ ] = { - 1 } ; int_T enableFcnCallFlag [
] = { 1 , 1 , 1 } ; const char * raccelLoadInputsAndAperiodicHitTimes (
SimStruct * S , const char * inportFileName , int * matFileFormat ) { return
rt_RAccelReadInportsMatFile ( S , inportFileName , matFileFormat ) ; }
#include "simstruc.h"
#include "fixedpoint.h"
#include "slsa_sim_engine.h"
#include "simtarget/slSimTgtSLExecSimBridge.h"
#define deuqmir1px (-1)
B rtB ; X rtX ; DW rtDW ; static SimStruct model_S ; SimStruct * const rtS =
& model_S ; void MdlInitialize ( void ) { rtX . o5hpeumlsb = rtP .
Integrator_IC ; rtX . jj3u5nl5xw = rtP . Integrator1_IC ; rtX . gjk4wr4fl2 =
rtP . Integrator_IC_i2wxdirj0b ; rtX . mmd3jbgtkz = rtP .
Integrator1_IC_oww15y3mqv ; rtX . lwhlamef2l = rtP .
Integrator1_IC_gd0uumlv2e ; rtDW . parg3eaod2 = ( rtInf ) ; rtDW . dejwxn1xy1
= ( rtInf ) ; rtX . logvqwtjje = rtP . Integrator_IC_gpeemvg0a0 ; rtX .
mwqg0fpixr = rtP . Integrator1_IC_epj0vehlxx ; rtX . fj4qnnj32q = rtP .
Integrator2_IC ; rtX . kup1qcjf04 = rtP . Integrator2_IC_cltf442w1a ; rtX .
asvy4gz4nf = rtP . Integrator2_IC_lurynlyjlu ; rtX . ncrlhqqtwv = rtP .
Integrator_IC_n5atnjvnfo ; rtX . b2dj4ozsrg = rtP . Integrator2_IC_dzcpbao03e
; rtDW . axuxe1t12v = deuqmir1px ; rtDW . ef11rieqzy = 0U ; rtDW . ekogk2wdou
= deuqmir1px ; rtDW . derxydwg12 = 0U ; rtDW . asp1sby2rd = deuqmir1px ; rtDW
. iqsjj14to5 = 0U ; rtDW . fpyg5ifl2j = deuqmir1px ; rtDW . aa1tc5fl4n = 0U ;
} void MdlStart ( void ) { { bool externalInputIsInDatasetFormat = false ;
void * pISigstreamManager = rt_GetISigstreamManager ( rtS ) ;
rtwISigstreamManagerGetInputIsInDatasetFormat ( pISigstreamManager , &
externalInputIsInDatasetFormat ) ; if ( externalInputIsInDatasetFormat ) { }
} { { { bool isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU
srcInfo ; sdiLabelU loggedName = sdiGetLabelFromChars ( "Modelo no lineal:8"
) ; sdiLabelU origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName
= sdiGetLabelFromChars ( "Modelo no lineal:8" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "sistema_completo/To Workspace" ) ; sdiLabelU blockSID
= sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( ""
) ; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars (
"Modelo no lineal:8" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. opv5nm4cbf . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"319d524c-8628-4f07-992a-d8b483e9062f" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . opv5nm4cbf . AQHandles , hDT , & srcInfo ) ; if ( rtDW . opv5nm4cbf .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . opv5nm4cbf . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
opv5nm4cbf . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . opv5nm4cbf .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . opv5nm4cbf . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . opv5nm4cbf . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . opv5nm4cbf . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"omega" ) ; sdiRegisterWksVariable ( rtDW . opv5nm4cbf . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Modelo no lineal:6" ) ;
sdiLabelU origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Modelo no lineal:6" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "sistema_completo/To Workspace1" ) ; sdiLabelU
blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Modelo no lineal:6" ) ; sdiAsyncRepoDataTypeHandle
hDT = sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; {
sdiComplexity sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. cqxqnkszbm . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"c9918524-2547-43fe-9c50-e8834edf33b3" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . cqxqnkszbm . AQHandles , hDT , & srcInfo ) ; if ( rtDW . cqxqnkszbm .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . cqxqnkszbm . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
cqxqnkszbm . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . cqxqnkszbm .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . cqxqnkszbm . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . cqxqnkszbm . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . cqxqnkszbm . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"T_m" ) ; sdiRegisterWksVariable ( rtDW . cqxqnkszbm . AQHandles , varName ,
"timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Mux" ) ; sdiLabelU origSigName
= sdiGetLabelFromChars ( "" ) ; sdiLabelU propName = sdiGetLabelFromChars (
"Mux" ) ; sdiLabelU blockPath = sdiGetLabelFromChars (
"sistema_completo/To Workspace2" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "Mux" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 3 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . dt1kwx4eay . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "e8e0ba06-057d-4cdd-b624-6df80b8cc842" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . dt1kwx4eay . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . dt1kwx4eay . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . dt1kwx4eay . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . dt1kwx4eay .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . dt1kwx4eay . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
dt1kwx4eay . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
dt1kwx4eay . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . dt1kwx4eay . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"i_abc" ) ; sdiRegisterWksVariable ( rtDW . dt1kwx4eay . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Mux1" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Mux1" ) ; sdiLabelU blockPath = sdiGetLabelFromChars
( "sistema_completo/To Workspace3" ) ; sdiLabelU blockSID =
sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath = sdiGetLabelFromChars ( "" )
; sdiDims sigDims ; sdiLabelU sigName = sdiGetLabelFromChars ( "Mux1" ) ;
sdiAsyncRepoDataTypeHandle hDT = sdiAsyncRepoGetBuiltInDataTypeHandle (
DATA_TYPE_DOUBLE ) ; { sdiComplexity sigComplexity = REAL ;
sdiSampleTimeContinuity stCont = SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray
[ 1 ] = { 3 } ; sigDims . nDims = 1 ; sigDims . dimensions = sigDimsArray ;
srcInfo . numBlockPathElems = 1 ; srcInfo . fullBlockPath = ( sdiFullBlkPathU
) & blockPath ; srcInfo . SID = ( sdiSignalIDU ) & blockSID ; srcInfo .
subPath = subPath ; srcInfo . portIndex = 0 + 1 ; srcInfo . signalName =
sigName ; srcInfo . sigSourceUUID = 0 ; rtDW . dnhkap2j1u . AQHandles =
sdiStartAsyncioQueueCreation ( hDT , & srcInfo , rt_dataMapInfo . mmi .
InstanceMap . fullPath , "c47c8c9b-54ae-4a60-8b90-9bbed4a84596" ,
sigComplexity , & sigDims , DIMENSIONS_MODE_FIXED , stCont , "" ) ;
sdiCompleteAsyncioQueueCreation ( rtDW . dnhkap2j1u . AQHandles , hDT , &
srcInfo ) ; if ( rtDW . dnhkap2j1u . AQHandles ) {
sdiSetSignalSampleTimeString ( rtDW . dnhkap2j1u . AQHandles , "Continuous" ,
0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW . dnhkap2j1u .
AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . dnhkap2j1u . AQHandles ,
ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings ( rtDW .
dnhkap2j1u . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName ( rtDW .
dnhkap2j1u . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . dnhkap2j1u . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"v_abc" ) ; sdiRegisterWksVariable ( rtDW . dnhkap2j1u . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "" ) ; sdiLabelU origSigName =
sdiGetLabelFromChars ( "" ) ; sdiLabelU propName = sdiGetLabelFromChars ( ""
) ; sdiLabelU blockPath = sdiGetLabelFromChars (
"sistema_completo/Observador de estados reducido/Gain3" ) ; sdiLabelU
blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 1 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. pqdji2nimh . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"826219a9-3961-46ff-80c9-8d8e68b7b32c" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . pqdji2nimh . AQHandles , hDT , & srcInfo ) ; if ( rtDW . pqdji2nimh .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . pqdji2nimh . AQHandles ,
"Continuous" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW .
pqdji2nimh . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . pqdji2nimh .
AQHandles , ssGetTaskTime ( rtS , 1 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . pqdji2nimh . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . pqdji2nimh . AQHandles , loggedName , origSigName , propName ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { } } } } { FWksInfo * fromwksInfo ; if ( (
fromwksInfo = ( FWksInfo * ) calloc ( 1 , sizeof ( FWksInfo ) ) ) == ( NULL )
) { ssSetErrorStatus ( rtS ,
"from workspace STRING(Name) memory allocation error" ) ; } else {
fromwksInfo -> origWorkspaceVarName =
 "Simulink.signaleditorblock.SimulationData.getData('c2lzdGVtYV9jb21wbGV0by9TaWduYWwgRWRpdG9yMg==','1')"
; fromwksInfo -> origDataTypeId = 0 ; fromwksInfo -> origIsComplex = 0 ;
fromwksInfo -> origWidth = 1 ; fromwksInfo -> origElSize = sizeof ( real_T )
; fromwksInfo -> data = ( void * ) rtP . fromWS_Signal1_Data0 ; fromwksInfo
-> nDataPoints = 2 ; fromwksInfo -> time = ( double * ) rtP .
fromWS_Signal1_Time0 ; rtDW . jrzquc22ff . TimePtr = fromwksInfo -> time ;
rtDW . jrzquc22ff . DataPtr = fromwksInfo -> data ; rtDW . jrzquc22ff .
RSimInfoPtr = fromwksInfo ; } rtDW . jptvtcgs5c . PrevIndex = 0 ; } {
FWksInfo * fromwksInfo ; if ( ( fromwksInfo = ( FWksInfo * ) calloc ( 1 ,
sizeof ( FWksInfo ) ) ) == ( NULL ) ) { ssSetErrorStatus ( rtS ,
"from workspace STRING(Name) memory allocation error" ) ; } else {
fromwksInfo -> origWorkspaceVarName =
 "Simulink.signaleditorblock.SimulationData.getData('c2lzdGVtYV9jb21wbGV0by9TaWduYWwgRWRpdG9yMw==','2')"
; fromwksInfo -> origDataTypeId = 0 ; fromwksInfo -> origIsComplex = 0 ;
fromwksInfo -> origWidth = 1 ; fromwksInfo -> origElSize = sizeof ( real_T )
; fromwksInfo -> data = ( void * ) rtP . FromWorkspace_Data0 ; fromwksInfo ->
nDataPoints = 141 ; fromwksInfo -> time = ( double * ) rtP .
FromWorkspace_Time0 ; rtDW . ml4gzok05f . TimePtr = fromwksInfo -> time ;
rtDW . ml4gzok05f . DataPtr = fromwksInfo -> data ; rtDW . ml4gzok05f .
RSimInfoPtr = fromwksInfo ; } rtDW . l2ca51dh03 . PrevIndex = 0 ; } {
FWksInfo * fromwksInfo ; if ( ( fromwksInfo = ( FWksInfo * ) calloc ( 1 ,
sizeof ( FWksInfo ) ) ) == ( NULL ) ) { ssSetErrorStatus ( rtS ,
"from workspace STRING(Name) memory allocation error" ) ; } else {
fromwksInfo -> origWorkspaceVarName =
 "Simulink.signaleditorblock.SimulationData.getData('c2lzdGVtYV9jb21wbGV0by9TaWduYWwgRWRpdG9y','2')"
; fromwksInfo -> origDataTypeId = 0 ; fromwksInfo -> origIsComplex = 0 ;
fromwksInfo -> origWidth = 1 ; fromwksInfo -> origElSize = sizeof ( real_T )
; fromwksInfo -> data = ( void * ) rtP . FromWorkspace_Data0_knbhfzr1k4 ;
fromwksInfo -> nDataPoints = 2 ; fromwksInfo -> time = ( double * ) rtP .
FromWorkspace_Time0_o2ayzzhu4m ; rtDW . diriss0ycm . TimePtr = fromwksInfo ->
time ; rtDW . diriss0ycm . DataPtr = fromwksInfo -> data ; rtDW . diriss0ycm
. RSimInfoPtr = fromwksInfo ; } rtDW . h33gajte4x . PrevIndex = 0 ; } {
FWksInfo * fromwksInfo ; if ( ( fromwksInfo = ( FWksInfo * ) calloc ( 1 ,
sizeof ( FWksInfo ) ) ) == ( NULL ) ) { ssSetErrorStatus ( rtS ,
"from workspace STRING(Name) memory allocation error" ) ; } else {
fromwksInfo -> origWorkspaceVarName =
 "Simulink.signaleditorblock.SimulationData.getData('c2lzdGVtYV9jb21wbGV0by9TaWduYWwgRWRpdG9yNA==','2')"
; fromwksInfo -> origDataTypeId = 0 ; fromwksInfo -> origIsComplex = 0 ;
fromwksInfo -> origWidth = 1 ; fromwksInfo -> origElSize = sizeof ( real_T )
; fromwksInfo -> data = ( void * ) rtP . FromWorkspace_Data0_ahxyd5aapb ;
fromwksInfo -> nDataPoints = 2 ; fromwksInfo -> time = ( double * ) rtP .
FromWorkspace_Time0_bh3fvyhkqh ; rtDW . ekcesp1k2p . TimePtr = fromwksInfo ->
time ; rtDW . ekcesp1k2p . DataPtr = fromwksInfo -> data ; rtDW . ekcesp1k2p
. RSimInfoPtr = fromwksInfo ; } rtDW . gjvhwykizl . PrevIndex = 0 ; } {
FWksInfo * fromwksInfo ; if ( ( fromwksInfo = ( FWksInfo * ) calloc ( 1 ,
sizeof ( FWksInfo ) ) ) == ( NULL ) ) { ssSetErrorStatus ( rtS ,
"from workspace STRING(Name) memory allocation error" ) ; } else {
fromwksInfo -> origWorkspaceVarName =
 "Simulink.signaleditorblock.SimulationData.getData('c2lzdGVtYV9jb21wbGV0by9UYW1i','2')"
; fromwksInfo -> origDataTypeId = 0 ; fromwksInfo -> origIsComplex = 0 ;
fromwksInfo -> origWidth = 1 ; fromwksInfo -> origElSize = sizeof ( real_T )
; fromwksInfo -> data = ( void * ) rtP . FromWorkspace_Data0_kevgbgk3mx ;
fromwksInfo -> nDataPoints = 2 ; fromwksInfo -> time = ( double * ) rtP .
FromWorkspace_Time0_fivuiqtoah ; rtDW . hrwve4sncv . TimePtr = fromwksInfo ->
time ; rtDW . hrwve4sncv . DataPtr = fromwksInfo -> data ; rtDW . hrwve4sncv
. RSimInfoPtr = fromwksInfo ; } rtDW . ojtfq2nxq4 . PrevIndex = 0 ; }
MdlInitialize ( ) ; } void MdlOutputs ( int_T tid ) { real_T ba1zh5u24z ;
real_T kqwxlrwn32 ; real_T pts2i1npil ; real_T nxw4oizztd ; real_T jh3h40fnu3
; real_T tmp [ 9 ] ; real_T abc [ 3 ] ; real_T lastTime ; real_T * lastU ;
int32_T i ; { real_T t = ssGetTaskTime ( rtS , 0 ) ; real_T * pTimeValues = (
real_T * ) rtDW . jrzquc22ff . TimePtr ; real_T * pDataValues = ( real_T * )
rtDW . jrzquc22ff . DataPtr ; int numPoints , lastPoint ; FWksInfo *
fromwksInfo = ( FWksInfo * ) rtDW . jrzquc22ff . RSimInfoPtr ; numPoints =
fromwksInfo -> nDataPoints ; lastPoint = numPoints - 1 ; if ( t < pTimeValues
[ 0 ] ) { ba1zh5u24z = 0.0 ; } else if ( t == pTimeValues [ lastPoint ] ) {
ba1zh5u24z = pDataValues [ lastPoint ] ; } else if ( t > pTimeValues [
lastPoint ] ) { ba1zh5u24z = 0.0 ; } else { int_T currTimeIndex = rtDW .
jptvtcgs5c . PrevIndex ; if ( t < pTimeValues [ currTimeIndex ] ) { while ( t
< pTimeValues [ currTimeIndex ] ) { currTimeIndex -- ; } } else { while ( t
>= pTimeValues [ currTimeIndex + 1 ] ) { currTimeIndex ++ ; } } ba1zh5u24z =
pDataValues [ currTimeIndex ] ; rtDW . jptvtcgs5c . PrevIndex = currTimeIndex
; } } rtB . k2sn1moatl = rtP . k_l / rtP . r * ba1zh5u24z ; rtB . m4zsvrd15a
= rtX . o5hpeumlsb ; rtB . jpgl5uqnen = 1.0 / rtP . r * rtB . m4zsvrd15a ;
rtB . hpaqwswri4 = rtB . k2sn1moatl * muDoubleScalarSin ( rtB . jpgl5uqnen )
; rtB . ojegsw1r4w = rtX . jj3u5nl5xw ; rtB . hpqeqr0o0i = rtP . K_sia * rtB
. ojegsw1r4w ; rtB . efmorfp5hk = rtX . gjk4wr4fl2 ; rtB . homtptr1wn = rtP .
K_sa * rtB . efmorfp5hk ; rtB . dbz2he40vm = rtX . mmd3jbgtkz ; rtB .
mbfoowmnkk = rtX . lwhlamef2l ; rtB . dlox32o0ua = rtP . Gain4_Gain * rtB .
mbfoowmnkk ; rtB . oh5qhuuk4l = rtB . dlox32o0ua - rtB . m4zsvrd15a ; rtB .
dx12nyer4s = rtP . k_theta2 * rtB . oh5qhuuk4l ; rtB . c5fdscj3uw = rtB .
dbz2he40vm + rtB . dx12nyer4s ; { real_T t = ssGetTaskTime ( rtS , 0 ) ;
real_T * pTimeValues = ( real_T * ) rtDW . ml4gzok05f . TimePtr ; real_T *
pDataValues = ( real_T * ) rtDW . ml4gzok05f . DataPtr ; int numPoints ,
lastPoint ; FWksInfo * fromwksInfo = ( FWksInfo * ) rtDW . ml4gzok05f .
RSimInfoPtr ; numPoints = fromwksInfo -> nDataPoints ; lastPoint = numPoints
- 1 ; if ( t < pTimeValues [ 0 ] ) { kqwxlrwn32 = 0.0 ; } else if ( t ==
pTimeValues [ lastPoint ] ) { kqwxlrwn32 = pDataValues [ lastPoint ] ; } else
if ( t > pTimeValues [ lastPoint ] ) { kqwxlrwn32 = 0.0 ; } else { int_T
currTimeIndex = rtDW . l2ca51dh03 . PrevIndex ; if ( t < pTimeValues [
currTimeIndex ] ) { while ( t < pTimeValues [ currTimeIndex ] ) {
currTimeIndex -- ; } } else { while ( t >= pTimeValues [ currTimeIndex + 1 ]
) { currTimeIndex ++ ; } } kqwxlrwn32 = pDataValues [ currTimeIndex ] ; rtDW
. l2ca51dh03 . PrevIndex = currTimeIndex ; } } rtB . lnhoszrvh2 = rtP . r *
kqwxlrwn32 ; { real_T t = ssGetTaskTime ( rtS , 0 ) ; real_T * pTimeValues =
( real_T * ) rtDW . diriss0ycm . TimePtr ; real_T * pDataValues = ( real_T *
) rtDW . diriss0ycm . DataPtr ; int numPoints , lastPoint ; FWksInfo *
fromwksInfo = ( FWksInfo * ) rtDW . diriss0ycm . RSimInfoPtr ; numPoints =
fromwksInfo -> nDataPoints ; lastPoint = numPoints - 1 ; if ( t < pTimeValues
[ 0 ] ) { pts2i1npil = 0.0 ; } else if ( t == pTimeValues [ lastPoint ] ) {
pts2i1npil = pDataValues [ lastPoint ] ; } else if ( t > pTimeValues [
lastPoint ] ) { pts2i1npil = 0.0 ; } else { int_T currTimeIndex = rtDW .
h33gajte4x . PrevIndex ; if ( t < pTimeValues [ currTimeIndex ] ) { while ( t
< pTimeValues [ currTimeIndex ] ) { currTimeIndex -- ; } } else { while ( t
>= pTimeValues [ currTimeIndex + 1 ] ) { currTimeIndex ++ ; } } pts2i1npil =
pDataValues [ currTimeIndex ] ; rtDW . h33gajte4x . PrevIndex = currTimeIndex
; } } rtB . c11k2clsrb = rtP . r * pts2i1npil ; if ( ( rtDW . parg3eaod2 >=
ssGetT ( rtS ) ) && ( rtDW . dejwxn1xy1 >= ssGetT ( rtS ) ) ) { rtB .
mqiqfyapeu = 0.0 ; } else { lastTime = rtDW . parg3eaod2 ; lastU = & rtDW .
na0wt4oepn ; if ( rtDW . parg3eaod2 < rtDW . dejwxn1xy1 ) { if ( rtDW .
dejwxn1xy1 < ssGetT ( rtS ) ) { lastTime = rtDW . dejwxn1xy1 ; lastU = & rtDW
. dcb22robm1 ; } } else if ( rtDW . parg3eaod2 >= ssGetT ( rtS ) ) { lastTime
= rtDW . dejwxn1xy1 ; lastU = & rtDW . dcb22robm1 ; } rtB . mqiqfyapeu = (
rtB . c11k2clsrb - * lastU ) / ( ssGetT ( rtS ) - lastTime ) ; } if ( rtP .
ManualSwitch1_CurrentSetting == 1 ) { rtB . kylrj4an0z = rtB . lnhoszrvh2 ; }
else { rtB . kylrj4an0z = rtB . mqiqfyapeu ; } rtB . dshebsf5x3 = rtB .
kylrj4an0z - rtB . c5fdscj3uw ; rtB . ly2xklr4fz = rtP . b_a * rtB .
dshebsf5x3 ; rtB . h1lgwz5ovu = ( rtB . hpqeqr0o0i + rtB . homtptr1wn ) + rtB
. ly2xklr4fz ; rtB . gifca5zh5d = rtP . b_eq * rtB . c5fdscj3uw ; rtB .
gysjfnsk1k = ( rtB . hpaqwswri4 + rtB . h1lgwz5ovu ) + rtB . gifca5zh5d ; rtB
. ajowt5z1q5 = rtX . logvqwtjje ; rtB . fwvyoitwh3 = rtX . mwqg0fpixr ; rtB .
netv40pb4f = rtX . fj4qnnj32q ; rtB . jo0ai5lhtc = rtP . Pp * rtB .
mbfoowmnkk ; rtDW . axuxe1t12v = deuqmir1px ; tmp [ 0 ] = muDoubleScalarCos (
rtB . jo0ai5lhtc ) ; tmp [ 3 ] = muDoubleScalarSin ( rtB . jo0ai5lhtc ) ; tmp
[ 6 ] = 1.0 ; tmp [ 1 ] = muDoubleScalarCos ( rtB . jo0ai5lhtc -
2.0943951023931953 ) ; tmp [ 4 ] = muDoubleScalarSin ( rtB . jo0ai5lhtc -
2.0943951023931953 ) ; tmp [ 7 ] = 1.0 ; tmp [ 2 ] = muDoubleScalarCos ( rtB
. jo0ai5lhtc + 2.0943951023931953 ) ; tmp [ 5 ] = muDoubleScalarSin ( rtB .
jo0ai5lhtc + 2.0943951023931953 ) ; tmp [ 8 ] = 1.0 ; for ( i = 0 ; i < 3 ; i
++ ) { abc [ i ] = ( tmp [ i + 3 ] * rtB . fwvyoitwh3 + tmp [ i ] * rtB .
ajowt5z1q5 ) + rtB . netv40pb4f ; } rtB . k532jfpdyn = rtP . Gain_Gain * abc
[ 0 ] ; rtB . nxola2u2g3 = rtP . Gain1_Gain * abc [ 1 ] ; rtB . kmls4v4rl0 =
rtP . Gain2_Gain * abc [ 2 ] ; rtB . af03zirre2 = rtP . Pp * rtB . dlox32o0ua
; rtDW . ekogk2wdou = deuqmir1px ; tmp [ 0 ] = 0.66666666666666663 *
muDoubleScalarCos ( rtB . af03zirre2 ) ; tmp [ 3 ] = muDoubleScalarCos ( rtB
. af03zirre2 - 2.0943951023931953 ) * 0.66666666666666663 ; tmp [ 6 ] =
muDoubleScalarCos ( rtB . af03zirre2 + 2.0943951023931953 ) *
0.66666666666666663 ; tmp [ 1 ] = 0.66666666666666663 * muDoubleScalarSin (
rtB . af03zirre2 ) ; tmp [ 4 ] = muDoubleScalarSin ( rtB . af03zirre2 -
2.0943951023931953 ) * 0.66666666666666663 ; tmp [ 7 ] = muDoubleScalarSin (
rtB . af03zirre2 + 2.0943951023931953 ) * 0.66666666666666663 ; tmp [ 2 ] =
0.33333333333333331 ; tmp [ 5 ] = 0.33333333333333331 ; tmp [ 8 ] =
0.33333333333333331 ; for ( i = 0 ; i < 3 ; i ++ ) { abc [ i ] = ( tmp [ i +
3 ] * rtB . nxola2u2g3 + tmp [ i ] * rtB . k532jfpdyn ) + tmp [ i + 6 ] * rtB
. kmls4v4rl0 ; } rtB . j0wuzgingg = abc [ 0 ] ; rtB . aduksi4izn = abc [ 1 ]
; rtB . pia0trm2jv = abc [ 2 ] ; rtB . dhufv0kct4 = ( rtP . L_d - rtP . L_q )
* rtB . aduksi4izn ; rtB . fd152jijlt = rtB . dhufv0kct4 + rtP . lambda_m ;
rtB . klsfhtexlj = 1.5 * rtP . Pp * rtB . fd152jijlt ; rtB . p3yvsgsymc = rtB
. gysjfnsk1k / rtB . klsfhtexlj ; rtB . elrtlardao = rtB . fwvyoitwh3 * rtB .
nf03pyptgh ; rtB . pxphckphjm = rtB . elrtlardao + rtP . lambda_m ; rtB .
k25kf23xp1 = rtB . ajowt5z1q5 * rtB . pxphckphjm ; rtB . g5whoifsex = 1.5 *
rtP . Pp * rtB . k25kf23xp1 ; rtB . fb0zj32vbu = rtX . kup1qcjf04 ; rtB .
kga0ttwxih = rtP . Gain3_Gain * rtB . fb0zj32vbu ; rtB . kf3sjuqaae = rtX .
asvy4gz4nf ; { if ( rtDW . opv5nm4cbf . AQHandles && ssGetLogOutput ( rtS ) )
{ sdiWriteSignal ( rtDW . opv5nm4cbf . AQHandles , ssGetTaskTime ( rtS , 0 )
, ( char * ) & rtB . kf3sjuqaae + 0 ) ; } } { if ( rtDW . cqxqnkszbm .
AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW . cqxqnkszbm .
AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB . g5whoifsex + 0 ) ;
} } rtB . fwd4ryyrfq [ 0 ] = rtB . k532jfpdyn ; rtB . fwd4ryyrfq [ 1 ] = rtB
. nxola2u2g3 ; rtB . fwd4ryyrfq [ 2 ] = rtB . kmls4v4rl0 ; { if ( rtDW .
dt1kwx4eay . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW .
dt1kwx4eay . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB .
fwd4ryyrfq [ 0 ] + 0 ) ; } } rtB . iamwwft4eh = rtB . p3yvsgsymc - rtB .
j0wuzgingg ; rtB . edszxpjalm = 5000.0 * rtP . L_q * rtB . iamwwft4eh ; rtB .
azzpx0o04y = rtP . Pp * rtB . c5fdscj3uw ; rtB . nd4l40ce0g = rtP . L_d * rtB
. aduksi4izn ; rtB . pr2t5fqsxs = rtB . nd4l40ce0g + rtP . lambda_m ; rtB .
nly0v00iut = rtB . azzpx0o04y * rtB . pr2t5fqsxs ; rtB . cri0rzhvjd = rtB .
fb0zj32vbu - rtP . T_s_ref ; rtB . o1eaam1gvh = rtP . alpha_Cu * rtB .
cri0rzhvjd ; rtB . i1cztljifv = rtB . o1eaam1gvh + rtP . Constant2_Value ;
rtB . hst2mc5opt = rtP . R_s_ref * rtB . i1cztljifv ; rtB . kpwiznui1d = rtB
. j0wuzgingg * rtB . hst2mc5opt ; rtB . ollu0ludgl = ( rtB . edszxpjalm + rtB
. nly0v00iut ) + rtB . kpwiznui1d ; rtB . bpf3qazr44 = rtB . kf3sjuqaae * rtB
. aduksi4izn ; rtB . nc1t1u3huh = rtP . L_d * rtP . Pp * rtB . bpf3qazr44 ;
rtB . dnzl4cfy5s = rtB . ollu0ludgl + rtB . nc1t1u3huh ; rtB . ougga4l4aa =
rtP . Constant_Value - rtB . aduksi4izn ; rtB . pxdilwla4w = 5000.0 * rtP .
L_d * rtB . ougga4l4aa ; rtB . gl1fdweeu4 = rtP . L_q * rtB . azzpx0o04y ;
rtB . fhjbzzgwbv = rtB . gl1fdweeu4 * rtB . j0wuzgingg ; rtB . kl3fuewomp =
rtB . aduksi4izn * rtB . hst2mc5opt ; rtB . h1yubcelv5 = ( rtB . pxdilwla4w -
rtB . fhjbzzgwbv ) + rtB . kl3fuewomp ; rtB . k1afv5mxmb = rtB . kf3sjuqaae *
rtB . j0wuzgingg ; rtB . lggf4lsbsd = - rtP . L_q * rtP . Pp * rtB .
k1afv5mxmb ; rtB . evhxefbhrl = rtB . h1yubcelv5 + rtB . lggf4lsbsd ; rtB .
gt3nq2cscb = rtB . pia0trm2jv * rtB . hst2mc5opt ; rtB . aqulfpaq14 = rtP .
Constant_Value - rtB . pia0trm2jv ; rtB . bd4jar1csl = 5000.0 * rtP . L_ls *
rtB . aqulfpaq14 ; rtB . dvu22fpphc = rtB . gt3nq2cscb + rtB . bd4jar1csl ;
rtDW . asp1sby2rd = deuqmir1px ; tmp [ 0 ] = muDoubleScalarCos ( rtB .
af03zirre2 ) ; tmp [ 3 ] = muDoubleScalarSin ( rtB . af03zirre2 ) ; tmp [ 6 ]
= 1.0 ; tmp [ 1 ] = muDoubleScalarCos ( rtB . af03zirre2 - 2.0943951023931953
) ; tmp [ 4 ] = muDoubleScalarSin ( rtB . af03zirre2 - 2.0943951023931953 ) ;
tmp [ 7 ] = 1.0 ; tmp [ 2 ] = muDoubleScalarCos ( rtB . af03zirre2 +
2.0943951023931953 ) ; tmp [ 5 ] = muDoubleScalarSin ( rtB . af03zirre2 +
2.0943951023931953 ) ; tmp [ 8 ] = 1.0 ; for ( i = 0 ; i < 3 ; i ++ ) { abc [
i ] = ( tmp [ i + 3 ] * rtB . evhxefbhrl + tmp [ i ] * rtB . dnzl4cfy5s ) +
rtB . dvu22fpphc ; } rtB . iss1dct3mr = abc [ 0 ] ; rtB . hpveetbx1j = abc [
1 ] ; rtB . eufxrbxe4z = abc [ 2 ] ; rtB . fnvmvvqodz [ 0 ] = rtB .
iss1dct3mr ; rtB . fnvmvvqodz [ 1 ] = rtB . hpveetbx1j ; rtB . fnvmvvqodz [ 2
] = rtB . eufxrbxe4z ; { if ( rtDW . dnhkap2j1u . AQHandles && ssGetLogOutput
( rtS ) ) { sdiWriteSignal ( rtDW . dnhkap2j1u . AQHandles , ssGetTaskTime (
rtS , 0 ) , ( char * ) & rtB . fnvmvvqodz [ 0 ] + 0 ) ; } } rtB . dz4hiszfdk
= rtX . ncrlhqqtwv ; if ( rtP . ManualSwitch_CurrentSetting == 1 ) { rtB .
lswxnkpgli = rtB . c11k2clsrb ; } else { rtB . lswxnkpgli = rtB . dz4hiszfdk
; } rtB . oqeh2ikwco = rtP . Pp * rtB . kf3sjuqaae ; rtB . jxvum3zvoh = rtP .
L_d * rtB . fwvyoitwh3 ; rtB . ou4cqm5an5 = rtP . lambda_m + rtB . jxvum3zvoh
; rtB . onoapnnemw = rtB . oqeh2ikwco * rtB . ou4cqm5an5 ; rtB . iktm1emn42 =
rtP . Gain_Gain_bpfw2sk42p * rtB . iss1dct3mr ; rtB . kwquioc0a2 = rtP .
Gain1_Gain_o051fegvds * rtB . hpveetbx1j ; rtB . caar3rhznz = rtP .
Gain2_Gain_fpelshvuq5 * rtB . eufxrbxe4z ; rtDW . fpyg5ifl2j = deuqmir1px ;
tmp [ 0 ] = 0.66666666666666663 * muDoubleScalarCos ( rtB . jo0ai5lhtc ) ;
tmp [ 3 ] = muDoubleScalarCos ( rtB . jo0ai5lhtc - 2.0943951023931953 ) *
0.66666666666666663 ; tmp [ 6 ] = muDoubleScalarCos ( rtB . jo0ai5lhtc +
2.0943951023931953 ) * 0.66666666666666663 ; tmp [ 1 ] = 0.66666666666666663
* muDoubleScalarSin ( rtB . jo0ai5lhtc ) ; tmp [ 4 ] = muDoubleScalarSin (
rtB . jo0ai5lhtc - 2.0943951023931953 ) * 0.66666666666666663 ; tmp [ 7 ] =
muDoubleScalarSin ( rtB . jo0ai5lhtc + 2.0943951023931953 ) *
0.66666666666666663 ; tmp [ 2 ] = 0.33333333333333331 ; tmp [ 5 ] =
0.33333333333333331 ; tmp [ 8 ] = 0.33333333333333331 ; for ( i = 0 ; i < 3 ;
i ++ ) { abc [ i ] = ( tmp [ i + 3 ] * rtB . kwquioc0a2 + tmp [ i ] * rtB .
iktm1emn42 ) + tmp [ i + 6 ] * rtB . caar3rhznz ; } rtB . gbixzpfrni = rtB .
hst2mc5opt * rtB . ajowt5z1q5 ; rtB . jwtzlohzhs = ( abc [ 0 ] - rtB .
onoapnnemw ) - rtB . gbixzpfrni ; rtB . cbv3kpprwx = 1.0 / rtP . L_q * rtB .
jwtzlohzhs ; rtB . cdjjap5tb4 = rtB . hst2mc5opt * rtB . fwvyoitwh3 ; rtB .
k5iqwa5dub = rtP . L_q * rtB . ajowt5z1q5 ; rtB . lbeji2oewc = rtB .
oqeh2ikwco * rtB . k5iqwa5dub ; rtB . crxm0t2udf = ( abc [ 1 ] - rtB .
cdjjap5tb4 ) + rtB . lbeji2oewc ; rtB . avres0c4zw = 1.0 / rtP . L_d * rtB .
crxm0t2udf ; rtB . buftnulajc = rtB . hst2mc5opt * rtB . netv40pb4f ; rtB .
pwtikpf0j0 = abc [ 2 ] - rtB . buftnulajc ; rtB . marzp3jauv = 1.0 / rtP .
L_ls * rtB . pwtikpf0j0 ; rtB . dymq05epmd = 1.0 / rtP . r * rtB . mbfoowmnkk
; rtB . d3woaso1vk = rtP . b_eq * rtB . kf3sjuqaae ; rtB . lata252i34 = rtP .
k_l * muDoubleScalarSin ( rtB . dymq05epmd ) ; { real_T t = ssGetTaskTime (
rtS , 0 ) ; real_T * pTimeValues = ( real_T * ) rtDW . ekcesp1k2p . TimePtr ;
real_T * pDataValues = ( real_T * ) rtDW . ekcesp1k2p . DataPtr ; int
numPoints , lastPoint ; FWksInfo * fromwksInfo = ( FWksInfo * ) rtDW .
ekcesp1k2p . RSimInfoPtr ; numPoints = fromwksInfo -> nDataPoints ; lastPoint
= numPoints - 1 ; if ( t < pTimeValues [ 0 ] ) { nxw4oizztd = 0.0 ; } else if
( t == pTimeValues [ lastPoint ] ) { nxw4oizztd = pDataValues [ lastPoint ] ;
} else if ( t > pTimeValues [ lastPoint ] ) { nxw4oizztd = 0.0 ; } else {
int_T currTimeIndex = rtDW . gjvhwykizl . PrevIndex ; if ( t < pTimeValues [
currTimeIndex ] ) { while ( t < pTimeValues [ currTimeIndex ] ) {
currTimeIndex -- ; } } else { while ( t >= pTimeValues [ currTimeIndex + 1 ]
) { currTimeIndex ++ ; } } nxw4oizztd = pDataValues [ currTimeIndex ] ; rtDW
. gjvhwykizl . PrevIndex = currTimeIndex ; } } rtB . p5w0nljbtt = ba1zh5u24z
* rtB . lata252i34 ; rtB . jkeaxbjqzp = nxw4oizztd + rtB . p5w0nljbtt ; rtB .
lpe0g0uuxy = 1.0 / rtP . r * rtB . jkeaxbjqzp ; rtB . kwzodnsn51 = ( rtB .
g5whoifsex - rtB . lpe0g0uuxy ) - rtB . d3woaso1vk ; rtB . kjg2pniu00 = 1.0 /
rtP . J_eq * rtB . kwzodnsn51 ; rtB . c4ogmaamk3 = rtB . netv40pb4f * rtB .
netv40pb4f ; rtB . jivucflgxz = rtP . Gain1_Gain_bfommv4lak * rtB .
c4ogmaamk3 ; rtB . b0djz0xuyx = rtP . Gain2_Gain_gt5xi45qct * rtB .
hst2mc5opt ; { real_T t = ssGetTaskTime ( rtS , 0 ) ; real_T * pTimeValues =
( real_T * ) rtDW . hrwve4sncv . TimePtr ; real_T * pDataValues = ( real_T *
) rtDW . hrwve4sncv . DataPtr ; int numPoints , lastPoint ; FWksInfo *
fromwksInfo = ( FWksInfo * ) rtDW . hrwve4sncv . RSimInfoPtr ; numPoints =
fromwksInfo -> nDataPoints ; lastPoint = numPoints - 1 ; if ( t < pTimeValues
[ 0 ] ) { jh3h40fnu3 = 0.0 ; } else if ( t == pTimeValues [ lastPoint ] ) {
jh3h40fnu3 = pDataValues [ lastPoint ] ; } else if ( t > pTimeValues [
lastPoint ] ) { jh3h40fnu3 = 0.0 ; } else { int_T currTimeIndex = rtDW .
ojtfq2nxq4 . PrevIndex ; if ( t < pTimeValues [ currTimeIndex ] ) { while ( t
< pTimeValues [ currTimeIndex ] ) { currTimeIndex -- ; } } else { while ( t
>= pTimeValues [ currTimeIndex + 1 ] ) { currTimeIndex ++ ; } } jh3h40fnu3 =
pDataValues [ currTimeIndex ] ; rtDW . ojtfq2nxq4 . PrevIndex = currTimeIndex
; } } rtB . moi1nwsely = rtB . fb0zj32vbu - jh3h40fnu3 ; rtB . apbcudoiru =
1.0 / rtP . R_ts * rtB . moi1nwsely ; rtB . bptosmp44x = rtB . ajowt5z1q5 *
rtB . ajowt5z1q5 ; rtB . fxltoioz2o = rtB . fwvyoitwh3 * rtB . fwvyoitwh3 ;
rtB . k3waxxng0u = ( rtB . bptosmp44x + rtB . fxltoioz2o ) + rtB . jivucflgxz
; rtB . ek2tbc2vjr = rtB . k3waxxng0u * rtB . b0djz0xuyx ; rtB . j1emwd2qae =
rtB . ek2tbc2vjr - rtB . apbcudoiru ; rtB . lsyqe5i55w = 1.0 / rtP . C_ts *
rtB . j1emwd2qae ; rtB . c4n4jb2dgo = rtX . b2dj4ozsrg ; rtB . ny2zw5nsgp =
rtP . k_i * rtB . c4n4jb2dgo ; if ( ssIsSampleHit ( rtS , 1 , 0 ) ) { { if (
rtDW . pqdji2nimh . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal (
rtDW . pqdji2nimh . AQHandles , ssGetTaskTime ( rtS , 1 ) , ( char * ) & rtB
. ny2zw5nsgp + 0 ) ; } } } rtB . fzxuttsywh = rtP . k_omega2 * rtB .
oh5qhuuk4l ; rtB . lb5025wpa3 = 1.0 / rtP . J_eq * rtB . h1lgwz5ovu ; rtB .
ib34wfcjm3 = ( rtB . ny2zw5nsgp + rtB . lb5025wpa3 ) + rtB . fzxuttsywh ;
UNUSED_PARAMETER ( tid ) ; } void MdlOutputsTID2 ( int_T tid ) { rtB .
nf03pyptgh = rtP . L_d - rtP . L_q ; UNUSED_PARAMETER ( tid ) ; } void
MdlUpdate ( int_T tid ) { real_T * lastU ; if ( rtDW . parg3eaod2 == ( rtInf
) ) { rtDW . parg3eaod2 = ssGetT ( rtS ) ; lastU = & rtDW . na0wt4oepn ; }
else if ( rtDW . dejwxn1xy1 == ( rtInf ) ) { rtDW . dejwxn1xy1 = ssGetT ( rtS
) ; lastU = & rtDW . dcb22robm1 ; } else if ( rtDW . parg3eaod2 < rtDW .
dejwxn1xy1 ) { rtDW . parg3eaod2 = ssGetT ( rtS ) ; lastU = & rtDW .
na0wt4oepn ; } else { rtDW . dejwxn1xy1 = ssGetT ( rtS ) ; lastU = & rtDW .
dcb22robm1 ; } * lastU = rtB . c11k2clsrb ; UNUSED_PARAMETER ( tid ) ; } void
MdlUpdateTID2 ( int_T tid ) { UNUSED_PARAMETER ( tid ) ; } void
MdlDerivatives ( void ) { XDot * _rtXdot ; _rtXdot = ( ( XDot * ) ssGetdX (
rtS ) ) ; _rtXdot -> o5hpeumlsb = rtB . c5fdscj3uw ; _rtXdot -> jj3u5nl5xw =
rtB . efmorfp5hk ; _rtXdot -> gjk4wr4fl2 = rtB . dshebsf5x3 ; _rtXdot ->
mmd3jbgtkz = rtB . ib34wfcjm3 ; _rtXdot -> lwhlamef2l = rtB . kf3sjuqaae ;
_rtXdot -> logvqwtjje = rtB . cbv3kpprwx ; _rtXdot -> mwqg0fpixr = rtB .
avres0c4zw ; _rtXdot -> fj4qnnj32q = rtB . marzp3jauv ; _rtXdot -> kup1qcjf04
= rtB . lsyqe5i55w ; _rtXdot -> asvy4gz4nf = rtB . kjg2pniu00 ; _rtXdot ->
ncrlhqqtwv = rtB . lnhoszrvh2 ; _rtXdot -> b2dj4ozsrg = rtB . oh5qhuuk4l ; }
void MdlProjection ( void ) { } void MdlTerminate ( void ) { rt_FREE ( rtDW .
jrzquc22ff . RSimInfoPtr ) ; rt_FREE ( rtDW . ml4gzok05f . RSimInfoPtr ) ;
rt_FREE ( rtDW . diriss0ycm . RSimInfoPtr ) ; rt_FREE ( rtDW . ekcesp1k2p .
RSimInfoPtr ) ; rt_FREE ( rtDW . hrwve4sncv . RSimInfoPtr ) ; { if ( rtDW .
opv5nm4cbf . AQHandles ) { sdiTerminateStreaming ( & rtDW . opv5nm4cbf .
AQHandles ) ; } } { if ( rtDW . cqxqnkszbm . AQHandles ) {
sdiTerminateStreaming ( & rtDW . cqxqnkszbm . AQHandles ) ; } } { if ( rtDW .
dt1kwx4eay . AQHandles ) { sdiTerminateStreaming ( & rtDW . dt1kwx4eay .
AQHandles ) ; } } { if ( rtDW . dnhkap2j1u . AQHandles ) {
sdiTerminateStreaming ( & rtDW . dnhkap2j1u . AQHandles ) ; } } { if ( rtDW .
pqdji2nimh . AQHandles ) { sdiTerminateStreaming ( & rtDW . pqdji2nimh .
AQHandles ) ; } } } static void mr_sistema_completo_cacheDataAsMxArray (
mxArray * destArray , mwIndex i , int j , const void * srcData , size_t
numBytes ) ; static void mr_sistema_completo_cacheDataAsMxArray ( mxArray *
destArray , mwIndex i , int j , const void * srcData , size_t numBytes ) {
mxArray * newArray = mxCreateUninitNumericMatrix ( ( size_t ) 1 , numBytes ,
mxUINT8_CLASS , mxREAL ) ; memcpy ( ( uint8_T * ) mxGetData ( newArray ) , (
const uint8_T * ) srcData , numBytes ) ; mxSetFieldByNumber ( destArray , i ,
j , newArray ) ; } static void mr_sistema_completo_restoreDataFromMxArray (
void * destData , const mxArray * srcArray , mwIndex i , int j , size_t
numBytes ) ; static void mr_sistema_completo_restoreDataFromMxArray ( void *
destData , const mxArray * srcArray , mwIndex i , int j , size_t numBytes ) {
memcpy ( ( uint8_T * ) destData , ( const uint8_T * ) mxGetData (
mxGetFieldByNumber ( srcArray , i , j ) ) , numBytes ) ; } static void
mr_sistema_completo_cacheBitFieldToMxArray ( mxArray * destArray , mwIndex i
, int j , uint_T bitVal ) ; static void
mr_sistema_completo_cacheBitFieldToMxArray ( mxArray * destArray , mwIndex i
, int j , uint_T bitVal ) { mxSetFieldByNumber ( destArray , i , j ,
mxCreateDoubleScalar ( ( real_T ) bitVal ) ) ; } static uint_T
mr_sistema_completo_extractBitFieldFromMxArray ( const mxArray * srcArray ,
mwIndex i , int j , uint_T numBits ) ; static uint_T
mr_sistema_completo_extractBitFieldFromMxArray ( const mxArray * srcArray ,
mwIndex i , int j , uint_T numBits ) { const uint_T varVal = ( uint_T )
mxGetScalar ( mxGetFieldByNumber ( srcArray , i , j ) ) ; return varVal & ( (
1u << numBits ) - 1u ) ; } static void
mr_sistema_completo_cacheDataToMxArrayWithOffset ( mxArray * destArray ,
mwIndex i , int j , mwIndex offset , const void * srcData , size_t numBytes )
; static void mr_sistema_completo_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) { uint8_T * varData = ( uint8_T * ) mxGetData (
mxGetFieldByNumber ( destArray , i , j ) ) ; memcpy ( ( uint8_T * ) & varData
[ offset * numBytes ] , ( const uint8_T * ) srcData , numBytes ) ; } static
void mr_sistema_completo_restoreDataFromMxArrayWithOffset ( void * destData ,
const mxArray * srcArray , mwIndex i , int j , mwIndex offset , size_t
numBytes ) ; static void mr_sistema_completo_restoreDataFromMxArrayWithOffset
( void * destData , const mxArray * srcArray , mwIndex i , int j , mwIndex
offset , size_t numBytes ) { const uint8_T * varData = ( const uint8_T * )
mxGetData ( mxGetFieldByNumber ( srcArray , i , j ) ) ; memcpy ( ( uint8_T *
) destData , ( const uint8_T * ) & varData [ offset * numBytes ] , numBytes )
; } static void mr_sistema_completo_cacheBitFieldToCellArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal )
; static void mr_sistema_completo_cacheBitFieldToCellArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal )
{ mxSetCell ( mxGetFieldByNumber ( destArray , i , j ) , offset ,
mxCreateDoubleScalar ( ( real_T ) fieldVal ) ) ; } static uint_T
mr_sistema_completo_extractBitFieldFromCellArrayWithOffset ( const mxArray *
srcArray , mwIndex i , int j , mwIndex offset , uint_T numBits ) ; static
uint_T mr_sistema_completo_extractBitFieldFromCellArrayWithOffset ( const
mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T numBits ) {
const uint_T fieldVal = ( uint_T ) mxGetScalar ( mxGetCell (
mxGetFieldByNumber ( srcArray , i , j ) , offset ) ) ; return fieldVal & ( (
1u << numBits ) - 1u ) ; } mxArray * mr_sistema_completo_GetDWork ( ) {
static const char_T * ssDWFieldNames [ 3 ] = { "rtB" , "rtDW" ,
"NULL_PrevZCX" , } ; mxArray * ssDW = mxCreateStructMatrix ( 1 , 1 , 3 ,
ssDWFieldNames ) ; mr_sistema_completo_cacheDataAsMxArray ( ssDW , 0 , 0 , (
const void * ) & ( rtB ) , sizeof ( rtB ) ) ; { static const char_T *
rtdwDataFieldNames [ 21 ] = { "rtDW.parg3eaod2" , "rtDW.na0wt4oepn" ,
"rtDW.dejwxn1xy1" , "rtDW.dcb22robm1" , "rtDW.asp1sby2rd" , "rtDW.ekogk2wdou"
, "rtDW.axuxe1t12v" , "rtDW.fpyg5ifl2j" , "rtDW.jptvtcgs5c" ,
"rtDW.l2ca51dh03" , "rtDW.h33gajte4x" , "rtDW.gjvhwykizl" , "rtDW.ojtfq2nxq4"
, "rtDW.iqsjj14to5" , "rtDW.derxydwg12" , "rtDW.ef11rieqzy" ,
"rtDW.aa1tc5fl4n" , "rtDW.hpxndlhxow" , "rtDW.d3jkgjc1oh" , "rtDW.gh0bm4uzee"
, "rtDW.lljtkeuzy4" , } ; mxArray * rtdwData = mxCreateStructMatrix ( 1 , 1 ,
21 , rtdwDataFieldNames ) ; mr_sistema_completo_cacheDataAsMxArray ( rtdwData
, 0 , 0 , ( const void * ) & ( rtDW . parg3eaod2 ) , sizeof ( rtDW .
parg3eaod2 ) ) ; mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 1 ,
( const void * ) & ( rtDW . na0wt4oepn ) , sizeof ( rtDW . na0wt4oepn ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 2 , ( const void * )
& ( rtDW . dejwxn1xy1 ) , sizeof ( rtDW . dejwxn1xy1 ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 3 , ( const void * )
& ( rtDW . dcb22robm1 ) , sizeof ( rtDW . dcb22robm1 ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 4 , ( const void * )
& ( rtDW . asp1sby2rd ) , sizeof ( rtDW . asp1sby2rd ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 5 , ( const void * )
& ( rtDW . ekogk2wdou ) , sizeof ( rtDW . ekogk2wdou ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 6 , ( const void * )
& ( rtDW . axuxe1t12v ) , sizeof ( rtDW . axuxe1t12v ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 7 , ( const void * )
& ( rtDW . fpyg5ifl2j ) , sizeof ( rtDW . fpyg5ifl2j ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 8 , ( const void * )
& ( rtDW . jptvtcgs5c ) , sizeof ( rtDW . jptvtcgs5c ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 9 , ( const void * )
& ( rtDW . l2ca51dh03 ) , sizeof ( rtDW . l2ca51dh03 ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 10 , ( const void * )
& ( rtDW . h33gajte4x ) , sizeof ( rtDW . h33gajte4x ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 11 , ( const void * )
& ( rtDW . gjvhwykizl ) , sizeof ( rtDW . gjvhwykizl ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 12 , ( const void * )
& ( rtDW . ojtfq2nxq4 ) , sizeof ( rtDW . ojtfq2nxq4 ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 13 , ( const void * )
& ( rtDW . iqsjj14to5 ) , sizeof ( rtDW . iqsjj14to5 ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 14 , ( const void * )
& ( rtDW . derxydwg12 ) , sizeof ( rtDW . derxydwg12 ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 15 , ( const void * )
& ( rtDW . ef11rieqzy ) , sizeof ( rtDW . ef11rieqzy ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 16 , ( const void * )
& ( rtDW . aa1tc5fl4n ) , sizeof ( rtDW . aa1tc5fl4n ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 17 , ( const void * )
& ( rtDW . hpxndlhxow ) , sizeof ( rtDW . hpxndlhxow ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 18 , ( const void * )
& ( rtDW . d3jkgjc1oh ) , sizeof ( rtDW . d3jkgjc1oh ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 19 , ( const void * )
& ( rtDW . gh0bm4uzee ) , sizeof ( rtDW . gh0bm4uzee ) ) ;
mr_sistema_completo_cacheDataAsMxArray ( rtdwData , 0 , 20 , ( const void * )
& ( rtDW . lljtkeuzy4 ) , sizeof ( rtDW . lljtkeuzy4 ) ) ; mxSetFieldByNumber
( ssDW , 0 , 1 , rtdwData ) ; } return ssDW ; } void
mr_sistema_completo_SetDWork ( const mxArray * ssDW ) { ( void ) ssDW ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtB ) , ssDW , 0
, 0 , sizeof ( rtB ) ) ; { const mxArray * rtdwData = mxGetFieldByNumber (
ssDW , 0 , 1 ) ; mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & (
rtDW . parg3eaod2 ) , rtdwData , 0 , 0 , sizeof ( rtDW . parg3eaod2 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . na0wt4oepn
) , rtdwData , 0 , 1 , sizeof ( rtDW . na0wt4oepn ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . dejwxn1xy1
) , rtdwData , 0 , 2 , sizeof ( rtDW . dejwxn1xy1 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . dcb22robm1
) , rtdwData , 0 , 3 , sizeof ( rtDW . dcb22robm1 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . asp1sby2rd
) , rtdwData , 0 , 4 , sizeof ( rtDW . asp1sby2rd ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . ekogk2wdou
) , rtdwData , 0 , 5 , sizeof ( rtDW . ekogk2wdou ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . axuxe1t12v
) , rtdwData , 0 , 6 , sizeof ( rtDW . axuxe1t12v ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . fpyg5ifl2j
) , rtdwData , 0 , 7 , sizeof ( rtDW . fpyg5ifl2j ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . jptvtcgs5c
) , rtdwData , 0 , 8 , sizeof ( rtDW . jptvtcgs5c ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . l2ca51dh03
) , rtdwData , 0 , 9 , sizeof ( rtDW . l2ca51dh03 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . h33gajte4x
) , rtdwData , 0 , 10 , sizeof ( rtDW . h33gajte4x ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . gjvhwykizl
) , rtdwData , 0 , 11 , sizeof ( rtDW . gjvhwykizl ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . ojtfq2nxq4
) , rtdwData , 0 , 12 , sizeof ( rtDW . ojtfq2nxq4 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . iqsjj14to5
) , rtdwData , 0 , 13 , sizeof ( rtDW . iqsjj14to5 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . derxydwg12
) , rtdwData , 0 , 14 , sizeof ( rtDW . derxydwg12 ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . ef11rieqzy
) , rtdwData , 0 , 15 , sizeof ( rtDW . ef11rieqzy ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . aa1tc5fl4n
) , rtdwData , 0 , 16 , sizeof ( rtDW . aa1tc5fl4n ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . hpxndlhxow
) , rtdwData , 0 , 17 , sizeof ( rtDW . hpxndlhxow ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . d3jkgjc1oh
) , rtdwData , 0 , 18 , sizeof ( rtDW . d3jkgjc1oh ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . gh0bm4uzee
) , rtdwData , 0 , 19 , sizeof ( rtDW . gh0bm4uzee ) ) ;
mr_sistema_completo_restoreDataFromMxArray ( ( void * ) & ( rtDW . lljtkeuzy4
) , rtdwData , 0 , 20 , sizeof ( rtDW . lljtkeuzy4 ) ) ; } } mxArray *
mr_sistema_completo_GetSimStateDisallowedBlocks ( ) { mxArray * data =
mxCreateCellMatrix ( 10 , 3 ) ; mwIndex subs [ 2 ] , offset ; { static const
char_T * blockType [ 10 ] = { "Scope" , "Scope" , "Scope" , "Scope" , "Scope"
, "Scope" , "Scope" , "Scope" , "Scope" , "Scope" , } ; static const char_T *
blockPath [ 10 ] = { "sistema_completo/Scope1" , "sistema_completo/Scope2" ,
"sistema_completo/Scope3" , "sistema_completo/T_m" , "sistema_completo/T_s" ,
"sistema_completo/i_abc" , "sistema_completo/i_qd0" ,
"sistema_completo/omega(m; m_obs; m*)" ,
"sistema_completo/theta(m; m_obs; m*)" , "sistema_completo/v_abc" , } ;
static const int reason [ 10 ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , }
; for ( subs [ 0 ] = 0 ; subs [ 0 ] < 10 ; ++ ( subs [ 0 ] ) ) { subs [ 1 ] =
0 ; offset = mxCalcSingleSubscript ( data , 2 , subs ) ; mxSetCell ( data ,
offset , mxCreateString ( blockType [ subs [ 0 ] ] ) ) ; subs [ 1 ] = 1 ;
offset = mxCalcSingleSubscript ( data , 2 , subs ) ; mxSetCell ( data ,
offset , mxCreateString ( blockPath [ subs [ 0 ] ] ) ) ; subs [ 1 ] = 2 ;
offset = mxCalcSingleSubscript ( data , 2 , subs ) ; mxSetCell ( data ,
offset , mxCreateDoubleScalar ( ( real_T ) reason [ subs [ 0 ] ] ) ) ; } }
return data ; } void MdlInitializeSizes ( void ) { ssSetNumContStates ( rtS ,
12 ) ; ssSetNumPeriodicContStates ( rtS , 0 ) ; ssSetNumY ( rtS , 0 ) ;
ssSetNumU ( rtS , 0 ) ; ssSetDirectFeedThrough ( rtS , 0 ) ;
ssSetNumSampleTimes ( rtS , 2 ) ; ssSetNumBlocks ( rtS , 156 ) ;
ssSetNumBlockIO ( rtS , 129 ) ; ssSetNumBlockParams ( rtS , 344 ) ; } void
MdlInitializeSampleTimes ( void ) { ssSetSampleTime ( rtS , 0 , 0.0 ) ;
ssSetSampleTime ( rtS , 1 , 0.0 ) ; ssSetOffsetTime ( rtS , 0 , 0.0 ) ;
ssSetOffsetTime ( rtS , 1 , 1.0 ) ; } void raccel_set_checksum ( ) {
ssSetChecksumVal ( rtS , 0 , 883753185U ) ; ssSetChecksumVal ( rtS , 1 ,
1892521921U ) ; ssSetChecksumVal ( rtS , 2 , 2740970809U ) ; ssSetChecksumVal
( rtS , 3 , 2593317770U ) ; }
#if defined(_MSC_VER)
#pragma optimize( "", off )
#endif
SimStruct * raccel_register_model ( ssExecutionInfo * executionInfo ) {
static struct _ssMdlInfo mdlInfo ; static struct _ssBlkInfo2 blkInfo2 ;
static struct _ssBlkInfoSLSize blkInfoSLSize ; ( void ) memset ( ( char_T * )
rtS , 0 , sizeof ( SimStruct ) ) ; ( void ) memset ( ( char_T * ) & mdlInfo ,
0 , sizeof ( struct _ssMdlInfo ) ) ; ( void ) memset ( ( char_T * ) &
blkInfo2 , 0 , sizeof ( struct _ssBlkInfo2 ) ) ; ( void ) memset ( ( char_T *
) & blkInfoSLSize , 0 , sizeof ( struct _ssBlkInfoSLSize ) ) ;
ssSetBlkInfo2Ptr ( rtS , & blkInfo2 ) ; ssSetBlkInfoSLSizePtr ( rtS , &
blkInfoSLSize ) ; ssSetMdlInfoPtr ( rtS , & mdlInfo ) ; ssSetExecutionInfo (
rtS , executionInfo ) ; slsaAllocOPModelData ( rtS ) ; { static time_T
mdlPeriod [ NSAMPLE_TIMES ] ; static time_T mdlOffset [ NSAMPLE_TIMES ] ;
static time_T mdlTaskTimes [ NSAMPLE_TIMES ] ; static int_T mdlTsMap [
NSAMPLE_TIMES ] ; static int_T mdlSampleHits [ NSAMPLE_TIMES ] ; static
boolean_T mdlTNextWasAdjustedPtr [ NSAMPLE_TIMES ] ; static int_T
mdlPerTaskSampleHits [ NSAMPLE_TIMES * NSAMPLE_TIMES ] ; static time_T
mdlTimeOfNextSampleHit [ NSAMPLE_TIMES ] ; { int_T i ; for ( i = 0 ; i <
NSAMPLE_TIMES ; i ++ ) { mdlPeriod [ i ] = 0.0 ; mdlOffset [ i ] = 0.0 ;
mdlTaskTimes [ i ] = 0.0 ; mdlTsMap [ i ] = i ; mdlSampleHits [ i ] = 1 ; } }
ssSetSampleTimePtr ( rtS , & mdlPeriod [ 0 ] ) ; ssSetOffsetTimePtr ( rtS , &
mdlOffset [ 0 ] ) ; ssSetSampleTimeTaskIDPtr ( rtS , & mdlTsMap [ 0 ] ) ;
ssSetTPtr ( rtS , & mdlTaskTimes [ 0 ] ) ; ssSetSampleHitPtr ( rtS , &
mdlSampleHits [ 0 ] ) ; ssSetTNextWasAdjustedPtr ( rtS , &
mdlTNextWasAdjustedPtr [ 0 ] ) ; ssSetPerTaskSampleHitsPtr ( rtS , &
mdlPerTaskSampleHits [ 0 ] ) ; ssSetTimeOfNextSampleHitPtr ( rtS , &
mdlTimeOfNextSampleHit [ 0 ] ) ; } ssSetSolverMode ( rtS ,
SOLVER_MODE_SINGLETASKING ) ; { ssSetBlockIO ( rtS , ( ( void * ) & rtB ) ) ;
( void ) memset ( ( ( void * ) & rtB ) , 0 , sizeof ( B ) ) ; } { real_T * x
= ( real_T * ) & rtX ; ssSetContStates ( rtS , x ) ; ( void ) memset ( ( void
* ) x , 0 , sizeof ( X ) ) ; } { void * dwork = ( void * ) & rtDW ;
ssSetRootDWork ( rtS , dwork ) ; ( void ) memset ( dwork , 0 , sizeof ( DW )
) ; } { static DataTypeTransInfo dtInfo ; ( void ) memset ( ( char_T * ) &
dtInfo , 0 , sizeof ( dtInfo ) ) ; ssSetModelMappingInfo ( rtS , & dtInfo ) ;
dtInfo . numDataTypes = 23 ; dtInfo . dataTypeSizes = & rtDataTypeSizes [ 0 ]
; dtInfo . dataTypeNames = & rtDataTypeNames [ 0 ] ; dtInfo . BTransTable = &
rtBTransTable ; dtInfo . PTransTable = & rtPTransTable ; dtInfo .
dataTypeInfoTable = rtDataTypeInfoTable ; }
sistema_completo_InitializeDataMapInfo ( ) ; ssSetIsRapidAcceleratorActive (
rtS , true ) ; ssSetRootSS ( rtS , rtS ) ; ssSetVersion ( rtS ,
SIMSTRUCT_VERSION_LEVEL2 ) ; ssSetModelName ( rtS , "sistema_completo" ) ;
ssSetPath ( rtS , "sistema_completo" ) ; ssSetTStart ( rtS , 0.0 ) ;
ssSetTFinal ( rtS , 300.0 ) ; { static RTWLogInfo rt_DataLoggingInfo ;
rt_DataLoggingInfo . loggingInterval = ( NULL ) ; ssSetRTWLogInfo ( rtS , &
rt_DataLoggingInfo ) ; } { { static int_T rt_LoggedStateWidths [ ] = { 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 } ; static int_T
rt_LoggedStateNumDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1
, 1 } ; static int_T rt_LoggedStateDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 ,
1 , 1 , 1 , 1 , 1 , 1 } ; static boolean_T rt_LoggedStateIsVarDims [ ] = { 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static BuiltInDTypeId
rt_LoggedStateDataTypeIds [ ] = { SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE } ; static int_T
rt_LoggedStateComplexSignals [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 } ; static RTWPreprocessingFcnPtr rt_LoggingStatePreprocessingFcnPtrs [
] = { ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , (
NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) } ; static
const char_T * rt_LoggedStateLabels [ ] = { "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" } ; static const char_T * rt_LoggedStateBlockNames [ ] =
{ "sistema_completo/Observador de estados reducido/Integrator" ,
"sistema_completo/Controlador PID/Integrator1" ,
"sistema_completo/Controlador PID/Integrator" ,
"sistema_completo/Observador de estados reducido/Integrator1" ,
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Integrator1" ,
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator" ,
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator1" ,
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator2" ,
"sistema_completo/Modelo no lineal/Subsistema_t&#xE9;rmico/Integrator2" ,
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Integrator2" ,
"sistema_completo/Set-point/Integrator" ,
"sistema_completo/Observador de estados reducido/Integrator2" } ; static
const char_T * rt_LoggedStateNames [ ] = { "" , "" , "" , "" , "" , "" , "" ,
"" , "" , "" , "" , "" } ; static boolean_T rt_LoggedStateCrossMdlRef [ ] = {
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static
RTWLogDataTypeConvert rt_RTWLogDataTypeConvert [ ] = { { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } } ; static int_T rt_LoggedStateIdxList [ ] = { 0 , 1 ,
2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 } ; static RTWLogSignalInfo
rt_LoggedStateSignalInfo = { 12 , rt_LoggedStateWidths ,
rt_LoggedStateNumDimensions , rt_LoggedStateDimensions ,
rt_LoggedStateIsVarDims , ( NULL ) , ( NULL ) , rt_LoggedStateDataTypeIds ,
rt_LoggedStateComplexSignals , ( NULL ) , rt_LoggingStatePreprocessingFcnPtrs
, { rt_LoggedStateLabels } , ( NULL ) , ( NULL ) , ( NULL ) , {
rt_LoggedStateBlockNames } , { rt_LoggedStateNames } ,
rt_LoggedStateCrossMdlRef , rt_RTWLogDataTypeConvert , rt_LoggedStateIdxList
} ; static void * rt_LoggedStateSignalPtrs [ 12 ] ; rtliSetLogXSignalPtrs (
ssGetRTWLogInfo ( rtS ) , ( LogSignalPtrsType ) rt_LoggedStateSignalPtrs ) ;
rtliSetLogXSignalInfo ( ssGetRTWLogInfo ( rtS ) , & rt_LoggedStateSignalInfo
) ; rt_LoggedStateSignalPtrs [ 0 ] = ( void * ) & rtX . o5hpeumlsb ;
rt_LoggedStateSignalPtrs [ 1 ] = ( void * ) & rtX . jj3u5nl5xw ;
rt_LoggedStateSignalPtrs [ 2 ] = ( void * ) & rtX . gjk4wr4fl2 ;
rt_LoggedStateSignalPtrs [ 3 ] = ( void * ) & rtX . mmd3jbgtkz ;
rt_LoggedStateSignalPtrs [ 4 ] = ( void * ) & rtX . lwhlamef2l ;
rt_LoggedStateSignalPtrs [ 5 ] = ( void * ) & rtX . logvqwtjje ;
rt_LoggedStateSignalPtrs [ 6 ] = ( void * ) & rtX . mwqg0fpixr ;
rt_LoggedStateSignalPtrs [ 7 ] = ( void * ) & rtX . fj4qnnj32q ;
rt_LoggedStateSignalPtrs [ 8 ] = ( void * ) & rtX . kup1qcjf04 ;
rt_LoggedStateSignalPtrs [ 9 ] = ( void * ) & rtX . asvy4gz4nf ;
rt_LoggedStateSignalPtrs [ 10 ] = ( void * ) & rtX . ncrlhqqtwv ;
rt_LoggedStateSignalPtrs [ 11 ] = ( void * ) & rtX . b2dj4ozsrg ; }
rtliSetLogT ( ssGetRTWLogInfo ( rtS ) , "tout" ) ; rtliSetLogX (
ssGetRTWLogInfo ( rtS ) , "" ) ; rtliSetLogXFinal ( ssGetRTWLogInfo ( rtS ) ,
"xFinal" ) ; rtliSetLogVarNameModifier ( ssGetRTWLogInfo ( rtS ) , "none" ) ;
rtliSetLogFormat ( ssGetRTWLogInfo ( rtS ) , 4 ) ; rtliSetLogMaxRows (
ssGetRTWLogInfo ( rtS ) , 0 ) ; rtliSetLogDecimation ( ssGetRTWLogInfo ( rtS
) , 1 ) ; rtliSetLogY ( ssGetRTWLogInfo ( rtS ) , "" ) ;
rtliSetLogYSignalInfo ( ssGetRTWLogInfo ( rtS ) , ( NULL ) ) ;
rtliSetLogYSignalPtrs ( ssGetRTWLogInfo ( rtS ) , ( NULL ) ) ; } { static
struct _ssStatesInfo2 statesInfo2 ; ssSetStatesInfo2 ( rtS , & statesInfo2 )
; } { static ssPeriodicStatesInfo periodicStatesInfo ;
ssSetPeriodicStatesInfo ( rtS , & periodicStatesInfo ) ; } { static
ssJacobianPerturbationBounds jacobianPerturbationBounds ;
ssSetJacobianPerturbationBounds ( rtS , & jacobianPerturbationBounds ) ; } {
static ssSolverInfo slvrInfo ; static boolean_T contStatesDisabled [ 12 ] ;
static real_T absTol [ 12 ] = { 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 ,
1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 } ; static
uint8_T absTolControl [ 12 ] = { 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U ,
0U , 0U , 0U } ; static real_T contStateJacPerturbBoundMinVec [ 12 ] ; static
real_T contStateJacPerturbBoundMaxVec [ 12 ] ; { int i ; for ( i = 0 ; i < 12
; ++ i ) { contStateJacPerturbBoundMinVec [ i ] = 0 ;
contStateJacPerturbBoundMaxVec [ i ] = rtGetInf ( ) ; } } ssSetSolverRelTol (
rtS , 0.001 ) ; ssSetStepSize ( rtS , 0.0 ) ; ssSetMinStepSize ( rtS , 0.0 )
; ssSetMaxNumMinSteps ( rtS , - 1 ) ; ssSetMinStepViolatedError ( rtS , 0 ) ;
ssSetMaxStepSize ( rtS , 0.001 ) ; ssSetSolverMaxOrder ( rtS , - 1 ) ;
ssSetSolverRefineFactor ( rtS , 1 ) ; ssSetOutputTimes ( rtS , ( NULL ) ) ;
ssSetNumOutputTimes ( rtS , 0 ) ; ssSetOutputTimesOnly ( rtS , 0 ) ;
ssSetOutputTimesIndex ( rtS , 0 ) ; ssSetZCCacheNeedsReset ( rtS , 0 ) ;
ssSetDerivCacheNeedsReset ( rtS , 0 ) ; ssSetNumNonContDerivSigInfos ( rtS ,
0 ) ; ssSetNonContDerivSigInfos ( rtS , ( NULL ) ) ; ssSetSolverInfo ( rtS ,
& slvrInfo ) ; ssSetSolverName ( rtS , "VariableStepAuto" ) ;
ssSetVariableStepSolver ( rtS , 1 ) ; ssSetSolverConsistencyChecking ( rtS ,
0 ) ; ssSetSolverAdaptiveZcDetection ( rtS , 0 ) ;
ssSetSolverRobustResetMethod ( rtS , 0 ) ; ssSetAbsTolVector ( rtS , absTol )
; ssSetAbsTolControlVector ( rtS , absTolControl ) ;
ssSetSolverAbsTol_Obsolete ( rtS , absTol ) ;
ssSetSolverAbsTolControl_Obsolete ( rtS , absTolControl ) ;
ssSetJacobianPerturbationBoundsMinVec ( rtS , contStateJacPerturbBoundMinVec
) ; ssSetJacobianPerturbationBoundsMaxVec ( rtS ,
contStateJacPerturbBoundMaxVec ) ; ssSetSolverStateProjection ( rtS , 0 ) ;
ssSetSolverMassMatrixType ( rtS , ( ssMatrixType ) 0 ) ;
ssSetSolverMassMatrixNzMax ( rtS , 0 ) ; ssSetModelOutputs ( rtS , MdlOutputs
) ; ssSetModelUpdate ( rtS , MdlUpdate ) ; ssSetModelDerivatives ( rtS ,
MdlDerivatives ) ; ssSetSolverMaxConsecutiveMinStep ( rtS , 1 ) ;
ssSetSolverShapePreserveControl ( rtS , 2 ) ; ssSetTNextTid ( rtS , INT_MIN )
; ssSetTNext ( rtS , rtMinusInf ) ; ssSetSolverNeedsReset ( rtS ) ;
ssSetNumNonsampledZCs ( rtS , 0 ) ; ssSetContStateDisabled ( rtS ,
contStatesDisabled ) ; ssSetSolverMaxConsecutiveMinStep ( rtS , 1 ) ; }
ssSetChecksumVal ( rtS , 0 , 883753185U ) ; ssSetChecksumVal ( rtS , 1 ,
1892521921U ) ; ssSetChecksumVal ( rtS , 2 , 2740970809U ) ; ssSetChecksumVal
( rtS , 3 , 2593317770U ) ; { static const sysRanDType rtAlwaysEnabled =
SUBSYS_RAN_BC_ENABLE ; static RTWExtModeInfo rt_ExtModeInfo ; static const
sysRanDType * systemRan [ 5 ] ; gblRTWExtModeInfo = & rt_ExtModeInfo ;
ssSetRTWExtModeInfo ( rtS , & rt_ExtModeInfo ) ;
rteiSetSubSystemActiveVectorAddresses ( & rt_ExtModeInfo , systemRan ) ;
systemRan [ 0 ] = & rtAlwaysEnabled ; systemRan [ 1 ] = & rtAlwaysEnabled ;
systemRan [ 2 ] = & rtAlwaysEnabled ; systemRan [ 3 ] = & rtAlwaysEnabled ;
systemRan [ 4 ] = & rtAlwaysEnabled ; rteiSetModelMappingInfoPtr (
ssGetRTWExtModeInfo ( rtS ) , & ssGetModelMappingInfo ( rtS ) ) ;
rteiSetChecksumsPtr ( ssGetRTWExtModeInfo ( rtS ) , ssGetChecksums ( rtS ) )
; rteiSetTPtr ( ssGetRTWExtModeInfo ( rtS ) , ssGetTPtr ( rtS ) ) ; }
slsaDisallowedBlocksForSimTargetOP ( rtS ,
mr_sistema_completo_GetSimStateDisallowedBlocks ) ;
slsaGetWorkFcnForSimTargetOP ( rtS , mr_sistema_completo_GetDWork ) ;
slsaSetWorkFcnForSimTargetOP ( rtS , mr_sistema_completo_SetDWork ) ;
rt_RapidReadMatFileAndUpdateParams ( rtS ) ; if ( ssGetErrorStatus ( rtS ) )
{ return rtS ; } return rtS ; }
#if defined(_MSC_VER)
#pragma optimize( "", on )
#endif
const int_T gblParameterTuningTid = 2 ; void MdlOutputsParameterSampleTime (
int_T tid ) { MdlOutputsTID2 ( tid ) ; }

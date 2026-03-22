#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "sistema_completo_capi_host.h"
#define sizeof(s) ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)
#ifndef SS_UINT64
#define SS_UINT64 17
#endif
#ifndef SS_INT64
#define SS_INT64 18
#endif
#else
#include "builtin_typeid_types.h"
#include "sistema_completo.h"
#include "sistema_completo_capi.h"
#include "sistema_completo_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s)               ((NULL))
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static const rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 3 , TARGET_STRING (
"sistema_completo/park_directa" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0
} , { 1 , 3 , TARGET_STRING ( "sistema_completo/park_directa" ) ,
TARGET_STRING ( "" ) , 1 , 0 , 0 , 0 , 0 } , { 2 , 3 , TARGET_STRING (
"sistema_completo/park_directa" ) , TARGET_STRING ( "" ) , 2 , 0 , 0 , 0 , 0
} , { 3 , 0 , TARGET_STRING (
"sistema_completo/park_directa/is_active_c2_sistema_completo" ) ,
TARGET_STRING ( "is_active_c2_sistema_completo" ) , 0 , 1 , 0 , 0 , 0 } , { 4
, 4 , TARGET_STRING ( "sistema_completo/park_inversa" ) , TARGET_STRING ( ""
) , 0 , 0 , 0 , 0 , 0 } , { 5 , 4 , TARGET_STRING (
"sistema_completo/park_inversa" ) , TARGET_STRING ( "" ) , 1 , 0 , 0 , 0 , 0
} , { 6 , 4 , TARGET_STRING ( "sistema_completo/park_inversa" ) ,
TARGET_STRING ( "" ) , 2 , 0 , 0 , 0 , 0 } , { 7 , 0 , TARGET_STRING (
"sistema_completo/park_inversa/is_active_c1_sistema_completo" ) ,
TARGET_STRING ( "is_active_c1_sistema_completo" ) , 0 , 1 , 0 , 0 , 0 } , { 8
, 0 , TARGET_STRING ( "sistema_completo/ " ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 9 , 0 , TARGET_STRING ( "sistema_completo/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 10 , 0 , TARGET_STRING (
"sistema_completo/Gain2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , {
11 , 0 , TARGET_STRING (
 "sistema_completo/TmpSignal ConversionAt_asyncqueue_inserted_for_To Workspace2Inport1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 12 , 0 , TARGET_STRING (
 "sistema_completo/TmpSignal ConversionAt_asyncqueue_inserted_for_To Workspace3Inport1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 13 , 0 , TARGET_STRING (
"sistema_completo/Product" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , {
14 , 0 , TARGET_STRING ( "sistema_completo/Product1" ) , TARGET_STRING ( "" )
, 0 , 0 , 0 , 0 , 0 } , { 15 , 0 , TARGET_STRING ( "sistema_completo/Sum1" )
, TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 16 , 0 , TARGET_STRING (
"sistema_completo/Sum3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 17
, 0 , TARGET_STRING ( "sistema_completo/Controlador PID/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 18 , 0 , TARGET_STRING (
"sistema_completo/Controlador PID/Gain2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 19 , 0 , TARGET_STRING (
"sistema_completo/Controlador PID/Gain3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 20 , 0 , TARGET_STRING (
"sistema_completo/Controlador PID/Integrator" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 21 , 0 , TARGET_STRING (
"sistema_completo/Controlador PID/Integrator1" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 22 , 0 , TARGET_STRING (
"sistema_completo/Controlador PID/Sum1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 23 , 0 , TARGET_STRING (
"sistema_completo/Controlador PID/Sum2" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 24 , 0 , TARGET_STRING ( "sistema_completo/Desacople/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 25 , 0 , TARGET_STRING (
"sistema_completo/Desacople/Gain3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 ,
0 } , { 26 , 0 , TARGET_STRING ( "sistema_completo/Desacople/Gain5" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 27 , 0 , TARGET_STRING (
"sistema_completo/Desacople/Product1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 ,
0 , 0 } , { 28 , 0 , TARGET_STRING ( "sistema_completo/Desacople/Product2" )
, TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 29 , 0 , TARGET_STRING (
"sistema_completo/Desacople/Product3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 ,
0 , 0 } , { 30 , 0 , TARGET_STRING ( "sistema_completo/Desacople/Product4" )
, TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 31 , 0 , TARGET_STRING (
"sistema_completo/Desacople/Product5" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 ,
0 , 0 } , { 32 , 0 , TARGET_STRING ( "sistema_completo/Desacople/Sum" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 33 , 0 , TARGET_STRING (
"sistema_completo/Desacople/Sum1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 ,
0 } , { 34 , 0 , TARGET_STRING ( "sistema_completo/Desacople/Sum2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 35 , 0 , TARGET_STRING (
"sistema_completo/Desacople/Sum3" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 ,
0 } , { 36 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/park_directa/is_active_c3_sistema_completo"
) , TARGET_STRING ( "is_active_c3_sistema_completo" ) , 0 , 1 , 0 , 0 , 0 } ,
{ 37 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/park_inversa/is_active_c4_sistema_completo"
) , TARGET_STRING ( "is_active_c4_sistema_completo" ) , 0 , 1 , 0 , 0 , 0 } ,
{ 38 , 0 , TARGET_STRING ( "sistema_completo/Modelo no lineal/Gain5" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 39 , 0 , TARGET_STRING (
"sistema_completo/Modulador de Tensión/Gain" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 40 , 0 , TARGET_STRING (
"sistema_completo/Modulador de Tensión/Gain1" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 41 , 0 , TARGET_STRING (
"sistema_completo/Modulador de Tensión/Gain2" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 42 , 0 , TARGET_STRING (
"sistema_completo/Modulador de corriente/Gain" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 43 , 0 , TARGET_STRING (
"sistema_completo/Modulador de corriente/Gain1" ) , TARGET_STRING ( "" ) , 0
, 0 , 0 , 0 , 0 } , { 44 , 0 , TARGET_STRING (
"sistema_completo/Modulador de corriente/Gain2" ) , TARGET_STRING ( "" ) , 0
, 0 , 0 , 0 , 0 } , { 45 , 0 , TARGET_STRING (
"sistema_completo/Modulador de corriente/Sum" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 46 , 0 , TARGET_STRING (
"sistema_completo/Modulador de corriente/Sum1" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 47 , 0 , TARGET_STRING (
"sistema_completo/Modulador de corriente/Sum2" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 48 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Gain1" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 49 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Gain2" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 50 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Gain3" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 51 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Gain4" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 52 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Gain5" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 53 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Divide" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 54 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Product" ) , TARGET_STRING ( "" ) , 0 ,
0 , 0 , 0 , 0 } , { 55 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Sum" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 56 , 0 , TARGET_STRING (
"sistema_completo/Modulador de torque/Sum1" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 57 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Gain" ) , TARGET_STRING ( ""
) , 0 , 0 , 0 , 0 , 0 } , { 58 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Gain1" ) , TARGET_STRING (
"" ) , 0 , 0 , 0 , 0 , 0 } , { 59 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Gain2" ) , TARGET_STRING (
"" ) , 0 , 0 , 0 , 0 , 0 } , { 60 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Gain3" ) , TARGET_STRING (
"" ) , 0 , 0 , 0 , 0 , 0 } , { 61 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Integrator" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 62 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Integrator1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 63 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Integrator2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 64 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Sum" ) , TARGET_STRING ( ""
) , 0 , 0 , 0 , 0 , 0 } , { 65 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Sum2" ) , TARGET_STRING ( ""
) , 0 , 0 , 0 , 0 , 0 } , { 66 , 0 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Sum3" ) , TARGET_STRING ( ""
) , 0 , 0 , 0 , 0 , 0 } , { 67 , 0 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 68 , 0 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain1" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 69 , 0 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain2" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 70 , 0 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain3" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 71 , 0 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain4" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 72 , 0 , TARGET_STRING (
"sistema_completo/Set-point/Derivative" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 73 , 0 , TARGET_STRING ( "sistema_completo/Set-point/Gain" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 74 , 0 , TARGET_STRING (
"sistema_completo/Set-point/Gain1" ) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 ,
0 } , { 75 , 0 , TARGET_STRING ( "sistema_completo/Set-point/Integrator" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 76 , 0 , TARGET_STRING (
"sistema_completo/Set-point/Manual Switch" ) , TARGET_STRING ( "" ) , 0 , 0 ,
0 , 0 , 0 } , { 77 , 0 , TARGET_STRING (
"sistema_completo/Set-point/Manual Switch1" ) , TARGET_STRING ( "" ) , 0 , 0
, 0 , 0 , 0 } , { 78 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain" ) ,
TARGET_STRING ( "i_qs'" ) , 0 , 0 , 0 , 0 , 0 } , { 79 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 80 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain2" ) ,
TARGET_STRING ( "i_ds'" ) , 0 , 0 , 0 , 0 , 0 } , { 81 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain3" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 82 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain4" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 83 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain5" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 84 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Gain6" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 85 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator" )
, TARGET_STRING ( "i_qs" ) , 0 , 0 , 0 , 0 , 0 } , { 86 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator1" )
, TARGET_STRING ( "i_ds" ) , 0 , 0 , 0 , 0 , 0 } , { 87 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator2" )
, TARGET_STRING ( "i_ds" ) , 0 , 0 , 0 , 0 , 0 } , { 88 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 89 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 90 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 91 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product3" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 92 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product4" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 93 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product5" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 94 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Product6" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 95 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Sum" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 96 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Sum1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 97 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Sum2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 98 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Sum3" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 99 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Sum4" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 100 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Sum5" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 1 } , { 101 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 102 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Gain2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 103 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Gain3" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 104 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Gain5" ) ,
TARGET_STRING ( "omega_m'" ) , 0 , 0 , 0 , 0 , 0 } , { 105 , 0 ,
TARGET_STRING ( "sistema_completo/Modelo no lineal/Subsistema_Mecanico/Gain6"
) , TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 106 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Integrator1" ) ,
TARGET_STRING ( "theta_m" ) , 0 , 0 , 0 , 0 , 0 } , { 107 , 0 , TARGET_STRING
( "sistema_completo/Modelo no lineal/Subsistema_Mecanico/Integrator2" ) ,
TARGET_STRING ( "omega_m" ) , 0 , 0 , 0 , 0 , 0 } , { 108 , 0 , TARGET_STRING
( "sistema_completo/Modelo no lineal/Subsistema_Mecanico/Product1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 109 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Sum1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 110 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Sum2" ) ,
TARGET_STRING ( "T_l" ) , 0 , 0 , 0 , 0 , 0 } , { 111 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 112 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 113 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain3" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 114 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain4" ) ,
TARGET_STRING ( "Rs" ) , 0 , 0 , 0 , 0 , 0 } , { 115 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain5" ) ,
TARGET_STRING ( "Ts'" ) , 0 , 0 , 0 , 0 , 0 } , { 116 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain6" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 117 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Integrator2" ) ,
TARGET_STRING ( "Ts" ) , 0 , 0 , 0 , 0 , 0 } , { 118 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Product1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 119 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Product2" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 120 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Product3" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 121 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Product4" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 0 , 0 , 0 } , { 122 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Sum1" ) , TARGET_STRING
( "" ) , 0 , 0 , 0 , 0 , 0 } , { 123 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Sum2" ) , TARGET_STRING
( "" ) , 0 , 0 , 0 , 0 , 0 } , { 124 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Sum3" ) , TARGET_STRING
( "" ) , 0 , 0 , 0 , 0 , 0 } , { 125 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Sum4" ) , TARGET_STRING
( "" ) , 0 , 0 , 0 , 0 , 0 } , { 126 , 0 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Sum5" ) , TARGET_STRING
( "" ) , 0 , 0 , 0 , 0 , 0 } , { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0
, 0 } } ; static const rtwCAPI_BlockParameters rtBlockParameters [ ] = { {
127 , TARGET_STRING ( "sistema_completo/Constant" ) , TARGET_STRING ( "Value"
) , 0 , 0 , 0 } , { 128 , TARGET_STRING (
"sistema_completo/Controlador PID/Integrator" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 129 , TARGET_STRING (
"sistema_completo/Controlador PID/Integrator1" ) , TARGET_STRING (
"InitialCondition" ) , 0 , 0 , 0 } , { 130 , TARGET_STRING (
"sistema_completo/Modulador de Tensión/Gain" ) , TARGET_STRING ( "Gain" ) , 0
, 0 , 0 } , { 131 , TARGET_STRING (
"sistema_completo/Modulador de Tensión/Gain1" ) , TARGET_STRING ( "Gain" ) ,
0 , 0 , 0 } , { 132 , TARGET_STRING (
"sistema_completo/Modulador de Tensión/Gain2" ) , TARGET_STRING ( "Gain" ) ,
0 , 0 , 0 } , { 133 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Integrator" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 134 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Integrator1" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 135 , TARGET_STRING (
"sistema_completo/Observador de estados reducido/Integrator2" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 136 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain" ) , TARGET_STRING ( "Gain" ) , 0 , 0
, 0 } , { 137 , TARGET_STRING ( "sistema_completo/Sensores ideales/Gain1" ) ,
TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 138 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain2" ) , TARGET_STRING ( "Gain" ) , 0 ,
0 , 0 } , { 139 , TARGET_STRING ( "sistema_completo/Sensores ideales/Gain3" )
, TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 140 , TARGET_STRING (
"sistema_completo/Sensores ideales/Gain4" ) , TARGET_STRING ( "Gain" ) , 0 ,
0 , 0 } , { 141 , TARGET_STRING ( "sistema_completo/Set-point/Integrator" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 142 , TARGET_STRING (
"sistema_completo/Set-point/Manual Switch" ) , TARGET_STRING (
"CurrentSetting" ) , 1 , 0 , 0 } , { 143 , TARGET_STRING (
"sistema_completo/Set-point/Manual Switch1" ) , TARGET_STRING (
"CurrentSetting" ) , 1 , 0 , 0 } , { 144 , TARGET_STRING (
"sistema_completo/Signal Editor/From Workspace" ) , TARGET_STRING ( "Time0" )
, 0 , 2 , 0 } , { 145 , TARGET_STRING (
"sistema_completo/Signal Editor/From Workspace" ) , TARGET_STRING ( "Data0" )
, 0 , 2 , 0 } , { 146 , TARGET_STRING (
"sistema_completo/Signal Editor2/fromWS_Signal 1" ) , TARGET_STRING ( "Time0"
) , 0 , 2 , 0 } , { 147 , TARGET_STRING (
"sistema_completo/Signal Editor2/fromWS_Signal 1" ) , TARGET_STRING ( "Data0"
) , 0 , 2 , 0 } , { 148 , TARGET_STRING (
"sistema_completo/Signal Editor3/From Workspace" ) , TARGET_STRING ( "Time0"
) , 0 , 3 , 0 } , { 149 , TARGET_STRING (
"sistema_completo/Signal Editor3/From Workspace" ) , TARGET_STRING ( "Data0"
) , 0 , 3 , 0 } , { 150 , TARGET_STRING (
"sistema_completo/Signal Editor4/From Workspace" ) , TARGET_STRING ( "Time0"
) , 0 , 2 , 0 } , { 151 , TARGET_STRING (
"sistema_completo/Signal Editor4/From Workspace" ) , TARGET_STRING ( "Data0"
) , 0 , 2 , 0 } , { 152 , TARGET_STRING (
"sistema_completo/Tamb/From Workspace" ) , TARGET_STRING ( "Time0" ) , 0 , 2
, 0 } , { 153 , TARGET_STRING ( "sistema_completo/Tamb/From Workspace" ) ,
TARGET_STRING ( "Data0" ) , 0 , 2 , 0 } , { 154 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator" )
, TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 155 , TARGET_STRING
( "sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator1"
) , TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 156 ,
TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Electromagnetico/Integrator2" )
, TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 157 , TARGET_STRING
( "sistema_completo/Modelo no lineal/Subsistema_Mecanico/Integrator1" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 158 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_Mecanico/Integrator2" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 159 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Constant2" ) ,
TARGET_STRING ( "Value" ) , 0 , 0 , 0 } , { 160 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain1" ) ,
TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 161 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Gain2" ) ,
TARGET_STRING ( "Gain" ) , 0 , 0 , 0 } , { 162 , TARGET_STRING (
"sistema_completo/Modelo no lineal/Subsistema_térmico/Integrator2" ) ,
TARGET_STRING ( "InitialCondition" ) , 0 , 0 , 0 } , { 0 , ( NULL ) , ( NULL
) , 0 , 0 , 0 } } ; static int_T rt_LoggedStateIdxList [ ] = { - 1 } ; static
const rtwCAPI_Signals rtRootInputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0
, 0 , 0 , 0 , 0 } } ; static const rtwCAPI_Signals rtRootOutputs [ ] = { { 0
, 0 , ( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const
rtwCAPI_ModelParameters rtModelParameters [ ] = { { 163 , TARGET_STRING (
"C_ts" ) , 0 , 0 , 0 } , { 164 , TARGET_STRING ( "J_eq" ) , 0 , 0 , 0 } , {
165 , TARGET_STRING ( "K_sa" ) , 0 , 0 , 0 } , { 166 , TARGET_STRING (
"K_sia" ) , 0 , 0 , 0 } , { 167 , TARGET_STRING ( "L_d" ) , 0 , 0 , 0 } , {
168 , TARGET_STRING ( "L_ls" ) , 0 , 0 , 0 } , { 169 , TARGET_STRING ( "L_q"
) , 0 , 0 , 0 } , { 170 , TARGET_STRING ( "Pp" ) , 0 , 0 , 0 } , { 171 ,
TARGET_STRING ( "R_s_ref" ) , 0 , 0 , 0 } , { 172 , TARGET_STRING ( "R_ts" )
, 0 , 0 , 0 } , { 173 , TARGET_STRING ( "T_s_ref" ) , 0 , 0 , 0 } , { 174 ,
TARGET_STRING ( "alpha_Cu" ) , 0 , 0 , 0 } , { 175 , TARGET_STRING ( "b_a" )
, 0 , 0 , 0 } , { 176 , TARGET_STRING ( "b_eq" ) , 0 , 0 , 0 } , { 177 ,
TARGET_STRING ( "k_i" ) , 0 , 0 , 0 } , { 178 , TARGET_STRING ( "k_l" ) , 0 ,
0 , 0 } , { 179 , TARGET_STRING ( "k_omega2" ) , 0 , 0 , 0 } , { 180 ,
TARGET_STRING ( "k_theta2" ) , 0 , 0 , 0 } , { 181 , TARGET_STRING (
"lambda_m" ) , 0 , 0 , 0 } , { 182 , TARGET_STRING ( "r" ) , 0 , 0 , 0 } , {
0 , ( NULL ) , 0 , 0 , 0 } } ;
#ifndef HOST_CAPI_BUILD
static void * rtDataAddrMap [ ] = { & rtB . j0wuzgingg , & rtB . aduksi4izn ,
& rtB . pia0trm2jv , & rtDW . derxydwg12 , & rtB . iss1dct3mr , & rtB .
hpveetbx1j , & rtB . eufxrbxe4z , & rtDW . iqsjj14to5 , & rtB . af03zirre2 ,
& rtB . lggf4lsbsd , & rtB . nc1t1u3huh , & rtB . fwd4ryyrfq [ 0 ] , & rtB .
fnvmvvqodz [ 0 ] , & rtB . k1afv5mxmb , & rtB . bpf3qazr44 , & rtB .
evhxefbhrl , & rtB . dnzl4cfy5s , & rtB . homtptr1wn , & rtB . hpqeqr0o0i , &
rtB . ly2xklr4fz , & rtB . efmorfp5hk , & rtB . ojegsw1r4w , & rtB .
dshebsf5x3 , & rtB . h1lgwz5ovu , & rtB . gl1fdweeu4 , & rtB . nd4l40ce0g , &
rtB . azzpx0o04y , & rtB . kpwiznui1d , & rtB . kl3fuewomp , & rtB .
fhjbzzgwbv , & rtB . nly0v00iut , & rtB . gt3nq2cscb , & rtB . pr2t5fqsxs , &
rtB . ollu0ludgl , & rtB . h1yubcelv5 , & rtB . dvu22fpphc , & rtDW .
aa1tc5fl4n , & rtDW . ef11rieqzy , & rtB . jo0ai5lhtc , & rtB . iktm1emn42 ,
& rtB . kwquioc0a2 , & rtB . caar3rhznz , & rtB . edszxpjalm , & rtB .
pxdilwla4w , & rtB . bd4jar1csl , & rtB . iamwwft4eh , & rtB . ougga4l4aa , &
rtB . aqulfpaq14 , & rtB . dhufv0kct4 , & rtB . gifca5zh5d , & rtB .
k2sn1moatl , & rtB . jpgl5uqnen , & rtB . klsfhtexlj , & rtB . p3yvsgsymc , &
rtB . hpaqwswri4 , & rtB . gysjfnsk1k , & rtB . fd152jijlt , & rtB .
fzxuttsywh , & rtB . dx12nyer4s , & rtB . lb5025wpa3 , & rtB . ny2zw5nsgp , &
rtB . m4zsvrd15a , & rtB . dbz2he40vm , & rtB . c4n4jb2dgo , & rtB .
oh5qhuuk4l , & rtB . ib34wfcjm3 , & rtB . c5fdscj3uw , & rtB . k532jfpdyn , &
rtB . nxola2u2g3 , & rtB . kmls4v4rl0 , & rtB . kga0ttwxih , & rtB .
dlox32o0ua , & rtB . mqiqfyapeu , & rtB . c11k2clsrb , & rtB . lnhoszrvh2 , &
rtB . dz4hiszfdk , & rtB . lswxnkpgli , & rtB . kylrj4an0z , & rtB .
cbv3kpprwx , & rtB . g5whoifsex , & rtB . avres0c4zw , & rtB . k5iqwa5dub , &
rtB . jxvum3zvoh , & rtB . oqeh2ikwco , & rtB . marzp3jauv , & rtB .
ajowt5z1q5 , & rtB . fwvyoitwh3 , & rtB . netv40pb4f , & rtB . k25kf23xp1 , &
rtB . cdjjap5tb4 , & rtB . lbeji2oewc , & rtB . gbixzpfrni , & rtB .
onoapnnemw , & rtB . buftnulajc , & rtB . elrtlardao , & rtB . crxm0t2udf , &
rtB . jwtzlohzhs , & rtB . ou4cqm5an5 , & rtB . pwtikpf0j0 , & rtB .
pxphckphjm , & rtB . nf03pyptgh , & rtB . d3woaso1vk , & rtB . lata252i34 , &
rtB . dymq05epmd , & rtB . kjg2pniu00 , & rtB . lpe0g0uuxy , & rtB .
mbfoowmnkk , & rtB . kf3sjuqaae , & rtB . p5w0nljbtt , & rtB . kwzodnsn51 , &
rtB . jkeaxbjqzp , & rtB . jivucflgxz , & rtB . b0djz0xuyx , & rtB .
o1eaam1gvh , & rtB . hst2mc5opt , & rtB . lsyqe5i55w , & rtB . apbcudoiru , &
rtB . fb0zj32vbu , & rtB . fxltoioz2o , & rtB . c4ogmaamk3 , & rtB .
bptosmp44x , & rtB . ek2tbc2vjr , & rtB . j1emwd2qae , & rtB . moi1nwsely , &
rtB . k3waxxng0u , & rtB . cri0rzhvjd , & rtB . i1cztljifv , & rtP .
Constant_Value , & rtP . Integrator_IC_i2wxdirj0b , & rtP . Integrator1_IC ,
& rtP . Gain_Gain_bpfw2sk42p , & rtP . Gain1_Gain_o051fegvds , & rtP .
Gain2_Gain_fpelshvuq5 , & rtP . Integrator_IC , & rtP .
Integrator1_IC_oww15y3mqv , & rtP . Integrator2_IC_dzcpbao03e , & rtP .
Gain_Gain , & rtP . Gain1_Gain , & rtP . Gain2_Gain , & rtP . Gain3_Gain , &
rtP . Gain4_Gain , & rtP . Integrator_IC_n5atnjvnfo , & rtP .
ManualSwitch_CurrentSetting , & rtP . ManualSwitch1_CurrentSetting , & rtP .
FromWorkspace_Time0_o2ayzzhu4m [ 0 ] , & rtP . FromWorkspace_Data0_knbhfzr1k4
[ 0 ] , & rtP . fromWS_Signal1_Time0 [ 0 ] , & rtP . fromWS_Signal1_Data0 [ 0
] , & rtP . FromWorkspace_Time0 [ 0 ] , & rtP . FromWorkspace_Data0 [ 0 ] , &
rtP . FromWorkspace_Time0_bh3fvyhkqh [ 0 ] , & rtP .
FromWorkspace_Data0_ahxyd5aapb [ 0 ] , & rtP . FromWorkspace_Time0_fivuiqtoah
[ 0 ] , & rtP . FromWorkspace_Data0_kevgbgk3mx [ 0 ] , & rtP .
Integrator_IC_gpeemvg0a0 , & rtP . Integrator1_IC_epj0vehlxx , & rtP .
Integrator2_IC , & rtP . Integrator1_IC_gd0uumlv2e , & rtP .
Integrator2_IC_lurynlyjlu , & rtP . Constant2_Value , & rtP .
Gain1_Gain_bfommv4lak , & rtP . Gain2_Gain_gt5xi45qct , & rtP .
Integrator2_IC_cltf442w1a , & rtP . C_ts , & rtP . J_eq , & rtP . K_sa , &
rtP . K_sia , & rtP . L_d , & rtP . L_ls , & rtP . L_q , & rtP . Pp , & rtP .
R_s_ref , & rtP . R_ts , & rtP . T_s_ref , & rtP . alpha_Cu , & rtP . b_a , &
rtP . b_eq , & rtP . k_i , & rtP . k_l , & rtP . k_omega2 , & rtP . k_theta2
, & rtP . lambda_m , & rtP . r , } ; static int32_T * rtVarDimsAddrMap [ ] =
{ ( NULL ) } ;
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , ( uint8_T ) SS_DOUBLE , 0 , 0 , 0 } ,
{ "unsigned char" , "uint8_T" , 0 , 0 , sizeof ( uint8_T ) , ( uint8_T )
SS_UINT8 , 0 , 0 , 0 } } ;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , } ; static const rtwCAPI_DimensionMap rtDimensionMap [ ] = { {
rtwCAPI_SCALAR , 0 , 2 , 0 } , { rtwCAPI_VECTOR , 2 , 2 , 0 } , {
rtwCAPI_VECTOR , 4 , 2 , 0 } , { rtwCAPI_VECTOR , 6 , 2 , 0 } } ; static
const uint_T rtDimensionArray [ ] = { 1 , 1 , 3 , 1 , 2 , 1 , 141 , 1 } ;
static const real_T rtcapiStoredFloats [ ] = { 0.0 } ; static const
rtwCAPI_FixPtMap rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) ,
rtwCAPI_FIX_RESERVED , 0 , 0 , ( boolean_T ) 0 } , } ; static const
rtwCAPI_SampleTimeMap rtSampleTimeMap [ ] = { { ( const void * ) &
rtcapiStoredFloats [ 0 ] , ( const void * ) & rtcapiStoredFloats [ 0 ] , (
int8_T ) 0 , ( uint8_T ) 0 } , { ( NULL ) , ( NULL ) , 2 , 0 } } ; static
rtwCAPI_ModelMappingStaticInfo mmiStatic = { { rtBlockSignals , 127 ,
rtRootInputs , 0 , rtRootOutputs , 0 } , { rtBlockParameters , 36 ,
rtModelParameters , 20 } , { ( NULL ) , 0 } , { rtDataTypeMap ,
rtDimensionMap , rtFixPtMap , rtElementMap , rtSampleTimeMap ,
rtDimensionArray } , "float" , { 883753185U , 1892521921U , 2740970809U ,
2593317770U } , ( NULL ) , 0 , ( boolean_T ) 0 , rt_LoggedStateIdxList } ;
const rtwCAPI_ModelMappingStaticInfo * sistema_completo_GetCAPIStaticMap (
void ) { return & mmiStatic ; }
#ifndef HOST_CAPI_BUILD
void sistema_completo_InitializeDataMapInfo ( void ) { rtwCAPI_SetVersion ( (
* rt_dataMapInfoPtr ) . mmi , 1 ) ; rtwCAPI_SetStaticMap ( ( *
rt_dataMapInfoPtr ) . mmi , & mmiStatic ) ; rtwCAPI_SetLoggingStaticMap ( ( *
rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ; rtwCAPI_SetDataAddressMap ( ( *
rt_dataMapInfoPtr ) . mmi , rtDataAddrMap ) ; rtwCAPI_SetVarDimsAddressMap (
( * rt_dataMapInfoPtr ) . mmi , rtVarDimsAddrMap ) ;
rtwCAPI_SetInstanceLoggingInfo ( ( * rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArray ( ( * rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( ( * rt_dataMapInfoPtr ) . mmi , 0 ) ; }
#else
#ifdef __cplusplus
extern "C" {
#endif
void sistema_completo_host_InitializeDataMapInfo (
sistema_completo_host_DataMapInfo_T * dataMap , const char * path ) {
rtwCAPI_SetVersion ( dataMap -> mmi , 1 ) ; rtwCAPI_SetStaticMap ( dataMap ->
mmi , & mmiStatic ) ; rtwCAPI_SetDataAddressMap ( dataMap -> mmi , ( NULL ) )
; rtwCAPI_SetVarDimsAddressMap ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetPath ( dataMap -> mmi , path ) ; rtwCAPI_SetFullPath ( dataMap ->
mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArray ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( dataMap -> mmi , 0 ) ; }
#ifdef __cplusplus
}
#endif
#endif

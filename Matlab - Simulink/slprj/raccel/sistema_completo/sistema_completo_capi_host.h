#ifndef RTW_HEADER_sistema_completo_cap_host_h__
#define RTW_HEADER_sistema_completo_cap_host_h__
#ifdef HOST_CAPI_BUILD
#include "rtw_capi.h"
#include "rtw_modelmap_simtarget.h"
typedef struct { rtwCAPI_ModelMappingInfo mmi ; }
sistema_completo_host_DataMapInfo_T ;
#ifdef __cplusplus
extern "C" {
#endif
void sistema_completo_host_InitializeDataMapInfo (
sistema_completo_host_DataMapInfo_T * dataMap , const char * path ) ;
#ifdef __cplusplus
}
#endif
#endif
#endif

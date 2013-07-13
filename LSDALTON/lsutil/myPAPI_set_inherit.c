// Set inherit structure for PAPI
// to ensure that FLOP counts are done correctly,
// also for OMP processes.
void mypapif_set_inherit_(int *es, int *retval)
{
#ifdef VAR_PAPI
#include <papi.h>
#include <string.h>

PAPI_option_t opt;

 *retval = PAPI_assign_eventset_component( *es, 0 );

 if (*retval != PAPI_OK) return;

 memset( &opt, 0x0, sizeof ( PAPI_option_t ) );
 opt.inherit.inherit = PAPI_INHERIT_ALL;
 opt.inherit.eventset = *es;
 
 *retval = PAPI_set_opt( PAPI_INHERIT, &opt ); 
#endif
}

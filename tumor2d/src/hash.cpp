#include <string.h>

unsigned long hash_float( float f )
{
    unsigned long ui;
    memcpy( &ui, &f, sizeof( float ) );
    return ui;
}

unsigned long hash_double( double f )
{
    unsigned long ui1, ui2;
    float *pf = (float*)&f;
    memcpy( &ui1, &pf[0], sizeof( float ) );
    memcpy( &ui2, &pf[1], sizeof( float ) );
    return ui1 + ui2;
}

unsigned long hash_array( double *v, int n){

	unsigned long hash_code = 0;
	for( int i=0; i<n; i++){
		//printf( "[ %e -> %u ]\n", v[i], hash_double( v[i]));
		hash_code = (hash_code + hash_double( v[i]));
	}

	return hash_code;
}

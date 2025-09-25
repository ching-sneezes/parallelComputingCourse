#ifndef WRITE_MATLAB_VECTOR_H
#define WRITE_MATLAB_VECTOR_H

// ---------------------------------------------------------------------------------------
// Save a vector to a matlab file (pointer version)
// ---------------------------------------------------------------------------------------
int writeMatlabVector( FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b )
{
    #define u(i) u_p[i-nd1a]
    const int numPerLine=8;
    fprintf(matlabFile,"%s=[",name);
    for( int i=nd1a; i<=nd1b; i++ )
    {
        fprintf(matlabFile,"%20.15e ",u(i));
        if( (i-nd1a) % numPerLine == numPerLine-1 )
            fprintf(matlabFile,"...\n");
    }
    fprintf(matlabFile,"];\n");
    return 0;
    #undef u
}

// ------------------------------------------------------------------------------------
// Save a vector to a matlab file (array version)
// ------------------------------------------------------------------------------------
int writeMatlabVector( FILE *matlabFile, RealArray & u, const char *name, int nd1a, int nd1b )
{
    const int numPerLine=8;
    fprintf(matlabFile,"%s=[",name);
    for( int i=nd1a; i<=nd1b; i++ )
    {
        fprintf(matlabFile,"%20.15e ",u(i));
        if( (i-nd1a) % numPerLine == numPerLine-1 )
            fprintf(matlabFile,"...\n");
    }
    fprintf(matlabFile,"];\n");
    return 0;
}

#endif

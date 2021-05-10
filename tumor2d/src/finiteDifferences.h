void initMatrix( VoronoiDiagram *voronoiDiagram);
void initMatrixImplicit( VoronoiDiagram *voronoiDiagram);
//double secondDerivationOxygen( VoronoiCell * cell);
double secondDerivationOxygen( VoronoiCell * cell, VoronoiDiagram *voronoiDiagram);
double UpdateSystem( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference);
double UpdateSystemImplicit( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference);
void ConjugateGradient( float **A, float *b, float *x, int N, int iterations);

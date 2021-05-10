%module nixTumor2d
%include <std_vector.i>
%include <typemaps.i>

%template(DoubleVector) std::vector<double>;

%{
extern void tumor2d_interface(double InitialRadius, double InitialQuiescentFraction, double MaxCellDivisionRate, double DivisionDepth,
						double ECMThresholdQuiescence, double ECMProductionRate, double ECMDegradationRate,
						double EndTime, double OutputRate, double profileTime, int profileDepth,int rand_seed,
						std::vector<double> &gc_out, std::vector<double> &ecm_out, std::vector<double> &prolif_out);
extern void tumor2d_default(std::vector<double> &gc_out, std::vector<double> &ecm_out, std::vector<double> &prolif_out);
%}

extern void tumor2d_interface(double InitialRadius, double InitialQuiescentFraction, double MaxCellDivisionRate, double DivisionDepth,
						double ECMThresholdQuiescence, double ECMProductionRate, double ECMDegradationRate,
						double EndTime, double OutputRate, double profileTime, int profileDepth,int rand_seed,
						std::vector<double> &OUTPUT, std::vector<double> &OUTPUT, std::vector<double> &OUTPUT);
extern void tumor2d_default(std::vector<double> &OUTPUT, std::vector<double> &OUTPUT, std::vector<double> &OUTPUT);
/* Copyright (c) 2003-2015 by Mike Jarvis
 *
 * TreeCorr is free software: redistribution and use in source and binary forms,
 * with or without modification, are permitted provided that the following
 * conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions, and the disclaimer given in the accompanying LICENSE
 *    file.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions, and the disclaimer given in the documentation
 *    and/or other materials provided with the distribution.
 */

#ifndef TreeCorr_Corr2_H
#define TreeCorr_Corr2_H

#include <vector>
#include <string>

#include "Cell.h"
#include "Field.h"

template <int DC1, int DC2>
struct XiData;

// BinnedCorr2 encapsulates a binned correlation function.
template <int DC1, int DC2>
class BinnedCorr2
{

public:

    BinnedCorr2(double minsep, double maxsep, int nbins, double binsize, double b,
                double* xi0, double* xi1, double* xi2, double* xi3,
                double* meanr, double* meanlogr, double* weight, double* npairs);
    BinnedCorr2(const BinnedCorr2& rhs, bool copy_data=true);
    ~BinnedCorr2();

    void clear();  // Set all data to 0.

    template <int M>
    void process(const Field<DC1, M>& field, bool dots);
    template <int M>
    void process(const Field<DC1, M>& field1, const Field<DC2, M>& field2, bool dots);
    template <int M>
    void processPairwise(const SimpleField<DC1, M>& field, const SimpleField<DC2, M>& field2,
                         bool dots);

    // Main worker functions for calculating the result
    template <int M>
    void process2(const Cell<DC1,M>& c12);

    template <int M>
    void process11(const Cell<DC1,M>& c1, const Cell<DC2,M>& c2);

    template <int M>
    void directProcess11(const Cell<DC1,M>& c1, const Cell<DC2,M>& c2, const double dsq);

    // Note: op= only copies _data.  Not all the params.
    void operator=(const BinnedCorr2<DC1,DC2>& rhs);
    void operator+=(const BinnedCorr2<DC1,DC2>& rhs);

protected:

    double _minsep;
    double _maxsep;
    int _nbins;
    double _binsize;
    double _b;
    double _logminsep;
    double _halfminsep;
    double _minsepsq;
    double _maxsepsq;
    double _bsq;
    int _metric; // Stores which Metric is being used for the analysis.

    // These are usually allocated in the python layer and just built up here.
    // So all we have here is a bare pointer for each of them.
    // However, for the OpenMP stuff, we do create copies that we need to delete.
    // So keep track if we own the data and need to delete the memory ourselves.
    bool _owns_data;

    // The different correlation functions have different numbers of arrays for xi, 
    // so encapsulate that difference with a templated XiData class.
    XiData<DC1,DC2> _xi;
    double* _meanr;
    double* _meanlogr;
    double* _weight;
    double* _npairs;
};

template <int DC1, int DC2>
struct XiData // This works for NK, KK
{
    XiData(double* xi0, double*, double*, double*) : xi(xi0) {}

    void new_data(int n) { xi = new double[n]; }
    void delete_data(int n) { delete [] xi; xi = 0; }
    void copy(const XiData<DC1,DC2>& rhs,int n) 
    { for (int i=0; i<n; ++i) xi[i] = rhs.xi[i]; }
    void add(const XiData<DC1,DC2>& rhs,int n) 
    { for (int i=0; i<n; ++i) xi[i] += rhs.xi[i]; }
    void clear(int n)
    { for (int i=0; i<n; ++i) xi[i] = 0.; }
    void write(std::ostream& os) const // Just used for debugging.  Print the first value.
    { os << xi[0]; }

    double* xi;
};

template <int DC1, int DC2>
inline std::ostream& operator<<(std::ostream& os, const XiData<DC1, DC2>& xi)
{ xi.write(os); return os; }

template <int DC1>
struct XiData<DC1, GData> // This works for NG, KG
{
    XiData(double* xi0, double* xi1, double*, double*) : xi(xi0), xi_im(xi1) {}

    void new_data(int n) 
    {
        xi = new double[n]; 
        xi_im = new double[n]; 
    }
    void delete_data(int n) 
    {
        delete [] xi; xi = 0; 
        delete [] xi_im; xi_im = 0; 
    }
    void copy(const XiData<DC1,GData>& rhs,int n) 
    { 
        for (int i=0; i<n; ++i) xi[i] = rhs.xi[i]; 
        for (int i=0; i<n; ++i) xi_im[i] = rhs.xi_im[i]; 
    }
    void add(const XiData<DC1,GData>& rhs,int n) 
    {
        for (int i=0; i<n; ++i) xi[i] += rhs.xi[i]; 
        for (int i=0; i<n; ++i) xi_im[i] += rhs.xi_im[i]; 
    }
    void clear(int n)
    { 
        for (int i=0; i<n; ++i) xi[i] = 0.;
        for (int i=0; i<n; ++i) xi_im[i] = 0.;
    }
    void write(std::ostream& os) const 
    { os << xi[0]<<','<<xi_im[0]; }

    double* xi;
    double* xi_im;
};

template <>
struct XiData<GData, GData>
{
    XiData(double* xi0, double* xi1, double* xi2, double* xi3) :
        xip(xi0), xip_im(xi1), xim(xi2), xim_im(xi3) {}

    void new_data(int n) 
    {
        xip = new double[n]; 
        xip_im = new double[n]; 
        xim = new double[n]; 
        xim_im = new double[n]; 
    }
    void delete_data(int n) 
    {
        delete [] xip; xip = 0; 
        delete [] xip_im; xip_im = 0; 
        delete [] xim; xim = 0; 
        delete [] xim_im; xim_im = 0; 
    }
    void copy(const XiData<GData,GData>& rhs,int n) 
    { 
        for (int i=0; i<n; ++i) xip[i] = rhs.xip[i]; 
        for (int i=0; i<n; ++i) xip_im[i] = rhs.xip_im[i]; 
        for (int i=0; i<n; ++i) xim[i] = rhs.xim[i]; 
        for (int i=0; i<n; ++i) xim_im[i] = rhs.xim_im[i]; 
    }
    void add(const XiData<GData,GData>& rhs,int n) 
    {
        for (int i=0; i<n; ++i) xip[i] += rhs.xip[i]; 
        for (int i=0; i<n; ++i) xip_im[i] += rhs.xip_im[i]; 
        for (int i=0; i<n; ++i) xim[i] += rhs.xim[i]; 
        for (int i=0; i<n; ++i) xim_im[i] += rhs.xim_im[i]; 
    }
    void clear(int n)
    { 
        for (int i=0; i<n; ++i) xip[i] = 0.;
        for (int i=0; i<n; ++i) xip_im[i] = 0.;
        for (int i=0; i<n; ++i) xim[i] = 0.;
        for (int i=0; i<n; ++i) xim_im[i] = 0.;
    }
    void write(std::ostream& os) const 
    { os << xip[0]<<','<<xip_im[0]<<','<<xim[0]<<','<<xim_im; }

    double* xip;
    double* xip_im;
    double* xim;
    double* xim_im;
};

template <>
struct XiData<NData, NData>
{
    XiData(double* , double* , double* , double* ) {}
    void new_data(int n) {}
    void delete_data(int n) {}
    void copy(const XiData<NData,NData>& rhs,int n) {}
    void add(const XiData<NData,NData>& rhs,int n) {}
    void clear(int n) {}
    void write(std::ostream& os) const {}
};


// The C interface for python
extern "C" {

    extern void* BuildNNCorr(double minsep, double maxsep, int nbins, double binsize, double b,
                             double* meanr, double* meanlogr, double* weight, double* npairs);
    extern void* BuildNKCorr(double minsep, double maxsep, int nbins, double binsize, double b,
                             double* xi,
                             double* meanr, double* meanlogr, double* weight, double* npairs);
    extern void* BuildNGCorr(double minsep, double maxsep, int nbins, double binsize, double b,
                             double* xi, double* xi_im,
                             double* meanr, double* meanlogr, double* weight, double* npairs);
    extern void* BuildKKCorr(double minsep, double maxsep, int nbins, double binsize, double b,
                             double* xi,
                             double* meanr, double* meanlogr, double* weight, double* npairs);
    extern void* BuildKGCorr(double minsep, double maxsep, int nbins, double binsize, double b,
                             double* xi, double* xi_im,
                             double* meanr, double* meanlogr, double* weight, double* npairs);
    extern void* BuildGGCorr(double minsep, double maxsep, int nbins, double binsize, double b,
                             double* xip, double* xip_im, double* xim, double* xim_im,
                             double* meanr, double* meanlogr, double* weight, double* npairs);

    extern void DestroyNNCorr(void* corr);
    extern void DestroyNKCorr(void* corr);
    extern void DestroyNGCorr(void* corr);
    extern void DestroyKKCorr(void* corr);
    extern void DestroyKGCorr(void* corr);
    extern void DestroyGGCorr(void* corr);

    extern void ProcessAutoNNFlat(void* corr, void* field, int dots);
    extern void ProcessAutoNN3D(void* corr, void* field, int dots);
    extern void ProcessAutoNNPerp(void* corr, void* field, int dots);
    extern void ProcessAutoKKFlat(void* corr, void* field, int dots);
    extern void ProcessAutoKK3D(void* corr, void* field, int dots);
    extern void ProcessAutoKKPerp(void* corr, void* field, int dots);
    extern void ProcessAutoGGFlat(void* corr, void* field, int dots);
    extern void ProcessAutoGG3D(void* corr, void* field, int dots);
    extern void ProcessAutoGGPerp(void* corr, void* field, int dots);

    extern void ProcessCrossNNFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNN3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNNPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNKFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNK3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNKPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNGFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNG3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossNGPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossKKFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossKK3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossKKPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossKGFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossKG3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossKGPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossGGFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossGG3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessCrossGGPerp(void* corr, void* field1, void* field2, int dots);

    extern void ProcessPairwiseNNFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNN3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNNPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNKFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNK3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNKPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNGFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNG3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseNGPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseKKFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseKK3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseKKPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseKGFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseKG3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseKGPerp(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseGGFlat(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseGG3D(void* corr, void* field1, void* field2, int dots);
    extern void ProcessPairwiseGGPerp(void* corr, void* field1, void* field2, int dots);

    extern int SetOMPThreads(int num_threads);
}

#endif

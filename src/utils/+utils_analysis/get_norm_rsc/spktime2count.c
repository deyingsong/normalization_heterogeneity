/*
 * SPKTIME2COUNT
 * Syntax (MATLAB):
 *   Y = spktime2count(s, neuronIdx, Tw, Ncount, option)
 * Description:
 *   Convert spike times to non-overlapping spike counts per neuron and time bin.
 *   - s         : 2 x Ns double matrix, columns [ts; id] with time in ms and integer neuron ID.
 *   - neuronIdx : vector of integer neuron IDs to include (order defines row order in Y).
 *   - Tw        : positive scalar, bin width in ms.
 *   - Ncount    : positive integer, number of bins (time from 0 to Tw*Ncount).
 *   - option    : 1 if neuronIdx is continuous & sorted (fast path), otherwise 0.
 * Outputs:
 *   - Y         : Nid x Ncount double matrix of spike counts.
 * Example:
 *   % s = [t; id] columns; ids 101..200; 50 ms bins, 100 bins
 *   Y = spktime2count(s, (101:200)', 50, 100, 1);
 *
 * Notes:
 *   - Robust to unsorted spike times; spikes outside [0, Tw*Ncount) are ignored.
 *   - Validates integer-valued neuron IDs.
 *   - For speed with non-contiguous IDs, provide sorted+contiguous IDs and set option=1.
 *   - Profile with: profile on; spktime2count(...); profile viewer
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <string.h>

/* -------- Utility: check near-integer -------- */
static int is_near_integer(double x){
    double r = floor(x + 0.5);
    return fabs(x - r) < 1e-9;
}

/* -------- Utility: is strictly increasing and contiguous -------- */
static int is_contiguous_sorted(const double *ids, mwSize n, int *first, int *last){
    mwSize i;
    if(n==0){ *first = 0; *last = -1; return 0; }
    if(!is_near_integer(ids[0])) return 0;
    int prev = (int)floor(ids[0] + 0.5);
    *first = prev;
    for(i=1;i<n;i++){
        if(!is_near_integer(ids[i])) return 0;
        int cur = (int)floor(ids[i] + 0.5);
        if(cur != prev + 1) return 0;
        prev = cur;
    }
    *last = prev;
    return 1;
}

/* -------- Utility: is strictly increasing (not necessarily contiguous) -------- */
static int is_strictly_increasing(const double *ids, mwSize n){
    mwSize i;
    if(n==0) return 0;
    int prev;
    if(!is_near_integer(ids[0])) return 0;
    prev = (int)floor(ids[0] + 0.5);
    for(i=1;i<n;i++){
        if(!is_near_integer(ids[i])) return 0;
        {
            int cur = (int)floor(ids[i] + 0.5);
            if(cur <= prev) return 0;
            prev = cur;
        }
    }
    return 1;
}

/* -------- Binary search on strictly increasing int array stored in double -------- */
static int bsearch_ids(const double *ids, mwSize n, int key){
    mwIndex lo = 0, hi = (mwIndex)n - 1;
    while(lo <= hi){
        mwIndex mid = lo + (hi - lo)/2;
        int v = (int)floor(ids[mid] + 0.5);
        if(v == key) return (int)mid;
        if(v < key) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

/* -------- Optional dense map for arbitrary IDs: maps [minID..maxID] -> row index or -1 -------- */
static int build_dense_map(const double *ids, mwSize n, int **map, int *minID, int *maxID){
    mwSize i; int mi, ma;
    if(n==0) return 0;
    mi = ma = (int)floor(ids[0] + 0.5);
    for(i=0;i<n;i++){
        if(!is_near_integer(ids[i])) return 0;
        {
            int v = (int)floor(ids[i] + 0.5);
            if(v < mi) mi = v; if(v > ma) ma = v;
        }
    }
    long range = (long)ma - (long)mi + 1L;
    if(range <= 0 || range > 10000000L || range > 20L*(long)n){
        /* Range too large for a dense map */
        return 0;
    }
    *map = (int*)mxCalloc((mwSize)range, sizeof(int));
    for(i=0;i<(mwSize)range;i++) (*map)[i] = -1;
    for(i=0;i<n;i++){
        int v = (int)floor(ids[i] + 0.5);
        (*map)[v - mi] = (int)i;
    }
    *minID = mi; *maxID = ma;
    return 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /* ---- Argument checks ---- */
    if(nrhs != 5){
        mexErrMsgIdAndTxt("spktime2count:arity", "Expected 5 inputs: s, neuronIdx, Tw, Ncount, option.");
    }
    if(nlhs > 1){
        mexErrMsgIdAndTxt("spktime2count:outputs", "One output (Y) expected.");
    }

    const mxArray *s_in = prhs[0];
    const mxArray *idx_in = prhs[1];

    if(!mxIsDouble(s_in) || mxIsComplex(s_in)){
        mexErrMsgIdAndTxt("spktime2count:sType", "s must be a real double matrix of size 2 x Ns.");
    }
    mwSize sM = mxGetM(s_in), sN = mxGetN(s_in);
    if(sM != 2){
        mexErrMsgIdAndTxt("spktime2count:sSize", "s must be 2 x Ns; got %u x %u.", (unsigned)sM, (unsigned)sN);
    }

    if(!mxIsDouble(idx_in) || mxIsComplex(idx_in)){
        mexErrMsgIdAndTxt("spktime2count:idxType", "neuronIdx must be a real double vector of integer IDs.");
    }
    mwSize idxM = mxGetM(idx_in), idxN = mxGetN(idx_in);
    if(idxM != 1 && idxN != 1){
        mexErrMsgIdAndTxt("spktime2count:idxSize", "neuronIdx must be a vector (Nx1 or 1xN).");
    }
    const double *ids = mxGetPr(idx_in);
    mwSize Nid = (idxM > idxN) ? idxM : idxN;

    double Tw = mxGetScalar(prhs[2]);
    if(!(Tw > 0 && mxIsDouble(prhs[2]))){
        mexErrMsgIdAndTxt("spktime2count:Tw", "Tw must be a positive real scalar.");
    }
    double Ncount_d = mxGetScalar(prhs[3]);
    if(!(Ncount_d > 0)){
        mexErrMsgIdAndTxt("spktime2count:Ncount", "Ncount must be a positive scalar integer.");
    }
    mwSize Ncount = (mwSize)floor(Ncount_d + 0.5);

    int option = (int)mxGetScalar(prhs[4]);

    /* ---- Precompute mapping strategy ---- */
    int ID1=0, ID2=-1, contiguous = 0, strict = 0, useDense = 0, minID=0, maxID=0;
    int *denseMap = NULL;
    contiguous = is_contiguous_sorted(ids, Nid, &ID1, &ID2);
    if(!contiguous){
        strict = is_strictly_increasing(ids, Nid);
        if(!strict){
            useDense = build_dense_map(ids, Nid, &denseMap, &minID, &maxID);
            if(!useDense){
                mexWarnMsgIdAndTxt("spktime2count:searchMode", "neuronIdx not sorted; falling back to linear search (may be slow). Consider sorting or using contiguous IDs with option=1.");
            }
        }
    }

    /* ---- Create output Y (Nid x Ncount) and zero it ---- */
    plhs[0] = mxCreateDoubleMatrix(Nid, Ncount, mxREAL);
    double *Y = mxGetPr(plhs[0]);
    memset(Y, 0, Nid * Ncount * sizeof(double));

    /* ---- Main accumulation ---- */
    const double *S = mxGetPr(s_in);
    mwSize Ns = sN;
    mwSize k;
    for(k=0; k<Ns; ++k){
        double ts = S[k*2 + 0];
        double idd = S[k*2 + 1];
        if(!(ts >= 0)) continue;               /* ignore negative times */
        mwSize bin = (mwSize)floor(ts / Tw);
        if(bin >= Ncount) continue;             /* ignore spikes beyond range */
        if(!is_near_integer(idd)) continue;     /* ignore non-integer IDs */
        int ID = (int)floor(idd + 0.5);

        /* Map neuron ID -> row j */
        int j = -1;
        if(option && contiguous){
            if(ID >= ID1 && ID <= ID2) j = ID - ID1;
        } else if(contiguous){
            if(ID >= ID1 && ID <= ID2) j = ID - ID1;
        } else if(strict){
            j = bsearch_ids(ids, Nid, ID);
        } else if(useDense){
            if(ID >= minID && ID <= maxID){
                j = denseMap[ID - minID];
            }
        } else {
            /* linear search */
            mwSize t;
            for(t=0; t<Nid; ++t){
                if((int)floor(ids[t] + 0.5) == ID){ j = (int)t; break; }
            }
        }
        if(j >= 0 && j < (int)Nid){
            Y[j + bin * Nid] += 1.0;
        }
    }

    if(denseMap) mxFree(denseMap);
}
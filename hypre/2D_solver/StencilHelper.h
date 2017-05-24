#ifndef STENCILHELPER_H
#define STENCILHELPER_H
class CoarseStencilHelper
{
	public:
	int currIndex[2];
	int nbrIndexCoarseLeft[2];
	int nbrIndexCoarseCenter[2];
	int nbrIndexCoarseRight[2];
	virtual ~CoarseStencilHelper(){};
	virtual int  nbr()                        = 0;
	virtual bool done()                       = 0;
	virtual CoarseStencilHelper &operator++() = 0;
};
class FineStencilHelper
{
	public:
	int currIndex[2];
	int nbrIndexFineLeftIn[2];
	int nbrIndexFineLeftOut[2];
	int nbrIndexFineRightIn[2];
	int nbrIndexFineRightOut[2];
	int nbrIndexFineOtherLeft[2];
	int nbrIndexFineOtherRight[2];
	virtual ~FineStencilHelper(){};
	virtual int  nbr()                      = 0;
	virtual bool done()                     = 0;
	virtual FineStencilHelper &operator++() = 0;
};

class EastFSH: public FineStencilHelper
{
    private:
	int     n;
	bool    left    = true;
	bool    is_done = false;
	Domain *d;

	public:
	EastFSH(Domain *d)
	{
		this->d                 = d;
		n                       = d->n;
		currIndex[0]            = n - 1;
		currIndex[1]            = n - 1;
		nbrIndexFineLeftIn[0]   = 0;
		nbrIndexFineLeftIn[1]   = n - 1;
		nbrIndexFineLeftOut[0]  = 1;
		nbrIndexFineLeftOut[1]  = n - 1;
		nbrIndexFineRightIn[0]  = 0;
		nbrIndexFineRightIn[1]  = n - 2;
		nbrIndexFineRightOut[0] = 1;
		nbrIndexFineRightOut[1] = n - 2;

		nbrIndexFineOtherLeft[0]  = n - 1;
		nbrIndexFineOtherLeft[1]  = n - 3;
		nbrIndexFineOtherRight[0] = n - 1;
		nbrIndexFineOtherRight[1] = 2;
	}
	int nbr()
	{
		if (left) {
			return d->nbr(Side::east);
		} else {
			return d->nbrRight(Side::east);
		}
	}
	EastFSH &operator++()
	{
		if (currIndex[1] > 0) {
			currIndex[1]--;
			nbrIndexFineLeftIn[1] -= 2;
			nbrIndexFineLeftOut[1] -= 2;
			nbrIndexFineRightIn[1] -= 2;
			nbrIndexFineRightOut[1] -= 2;
			if (nbrIndexFineLeftIn[1] < 0) {
				left = false;
				nbrIndexFineLeftIn[1] += n;
				nbrIndexFineLeftOut[1] += n;
				nbrIndexFineRightIn[1] += n;
				nbrIndexFineRightOut[1] += n;
			}
		}else{
            is_done=true;
        }
		return *this;
	}
    bool done(){
        return is_done;
    }
};
class WestCSH: public CoarseStencilHelper
{
    private:
	int     n;
	bool    first = false;
	bool    last  = false;
	bool    left  = true;
	bool    is_done = false;
	bool    coarse_left;
	Domain *d;
	public:
	WestCSH(Domain* d){
		this->d     = d;
		n           = d->n;
		coarse_left = d->isCoarseLeft(Side::west);
		if (coarse_left) {
			nbrIndexCoarseLeft[0]   = n - 1;
			nbrIndexCoarseLeft[1]   = n / 2 - 1;
			nbrIndexCoarseCenter[0] = n - 1;
			nbrIndexCoarseCenter[1] = n / 2;
			nbrIndexCoarseRight[0]  = n - 1;
			nbrIndexCoarseRight[1]  = n / 2 + 1;
		} else {
			first                   = true;
			nbrIndexCoarseLeft[0]   = n - 1;
			nbrIndexCoarseLeft[1]   = 0;
			nbrIndexCoarseCenter[0] = n - 1;
			nbrIndexCoarseCenter[1] = 1;
			nbrIndexCoarseRight[0]  = n - 1;
			nbrIndexCoarseRight[1]  = 2;
		}

		currIndex[0] = 0;
		currIndex[1] = 0;
	}
	int  nbr() { return d->nbr(Side::west); }
	bool done() { return is_done; }
	WestCSH &operator++()
	{
		if (currIndex[1] < n - 1) {
			if (first) {
				if (left) {
					left = false;
				} else {
					first = false;
					left  = true;
				}
			} else if (!last) {
				if (left) {
					left = false;
				} else {
					left = true;
					nbrIndexCoarseLeft[1]++;
					nbrIndexCoarseCenter[1]++;
					nbrIndexCoarseRight[1]++;
				}
			}
			currIndex[1]++;
			if (coarse_left && (currIndex[1] == n - 3)) {
				last = true;
			}
		} else {
			is_done = true;
		}
		return *this;
	}
};

/*
class NorthHelper: public StencilHelper
{
    private:
    std::array<int, 2> curr_index;
    int     i     = 0;
    int     n;
    bool    first = true;
    bool    last  = false;
    Domain *d;
    bool    left() { return curr_index[0] < d->n / 2; }
    public:
    NorthHelper(Domain* d){
        this->d       = d;
        n             = d->n;
        curr_index[0] = 0;
        curr_index[0] = n - 1;
    }
    std::array<int, 2> currIndex() {return curr_index;};
    int nbr()
    {
        if (left()) {
            return d->nbr(Side::north);
        } else {
            return d->nbrRight(Side::north);
        }
    }
    std::array<int, 2> nbrIndexFineLeftIn(){
        if (left()) {
            std::array<int, 2> ret = {n - 2*i, 0};
            return ret;
        } else {
            std::array<int, 2> ret = {2 * n - 2 * i, 0};
            return ret;
        }
    }
    std::array<int, 2> nbrIndexFineLeftOut(){
        if (left()) {
            std::array<int, 2> ret = {n - 2*i, n - 2};
            return ret;
        } else {
            std::array<int, 2> ret = {2 * n - 2 * i, n - 2};
            return ret;
        }
    }
    std::array<int, 2> nbrIndexFineRightIn(){
        if (left()) {
            std::array<int, 2> ret = {n - 2 * i - 1, n - 2};
            return ret;
        } else {
            std::array<int, 2> ret = {2 * n - 2 * i - 1, n - 2};
            return ret;
        }
    }
    std::array<int, 2> nbrIndexFineRightOut()   = 0;
    NorthHelper &operator++()
    {
        if (curr_index < n - 1) {
            curr_index[0]++;
            i++;
        }
        return *this;
    }
};
*/
#endif

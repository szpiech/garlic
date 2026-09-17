#ifndef __GARLIC_MATRIX_H__
#define __GARLIC_MATRIX_H__

#include <vector>
#include <cstring>

//A row-indexed matrix held in ONE contiguous buffer.
//
//This replaces the `T **rows` + `new T[nrow*ncol]` + row-pointer fixup pattern
//that initHapData/initWinData/initGLData/initLDData all hand-rolled, and the
//paired releaseX functions that had to remember to free row 0 and the pointer
//array but not each row.  Contiguity is deliberate and must be preserved: it
//is what made the LOD stage cache-friendly, and dropping it would undo that.
//
//operator[] returns T*, so existing `m[i][j]` call sites do not change -- that
//is the whole reason this is a header change rather than a sweep of every use.
//
//NOTE ON bool: do NOT instantiate this as Matrix<bool>.  std::vector<bool> is
//the proxy specialisation, has no T* to hand out, and operator[] would not
//compile; "fixing" that by returning a proxy would silently change what every
//firstCopy[i][j] site means.  Use unsigned char.
template <typename T>
class Matrix
{
public:
    Matrix() : nrow_(0), ncol_(0) {}

    //Copy and assignment MUST rebuild the row index.  The compiler-generated
    //versions would copy row pointers that aim into the SOURCE's buffer, which
    //is a dangling-pointer bug of exactly the shape this refactor exists to
    //remove (see B1).  Hence these are written out rather than defaulted.
    Matrix(const Matrix &o) : buf_(o.buf_), nrow_(o.nrow_), ncol_(o.ncol_) { repoint(); }
    Matrix &operator=(const Matrix &o)
    {
        if (this != &o) { buf_ = o.buf_; nrow_ = o.nrow_; ncol_ = o.ncol_; repoint(); }
        return *this;
    }

    void resize(long nrow, long ncol)
    {
        buf_.resize(size_t(nrow) * size_t(ncol));
        nrow_ = nrow; ncol_ = ncol;
        repoint();
    }

    //Sized and value-initialised in one call.
    void assign(long nrow, long ncol, const T &fill)
    {
        buf_.assign(size_t(nrow) * size_t(ncol), fill);
        nrow_ = nrow; ncol_ = ncol;
        repoint();
    }

    void clear() { buf_.clear(); row_.clear(); nrow_ = 0; ncol_ = 0; }

    //Row-at-a-time construction, for the reader path that gathers per-locus
    //rows it allocated before it knew the locus count.  resize()/assign() would
    //zero-fill the whole block first, committing every page up front; the old
    //new[] committed them progressively as the memcpy touched them, which
    //overlapped with freeing the reader's rows and so had a lower peak.  This
    //reproduces that: reserve once, then copy each row in with no zero pass.
    void reserveRows(long nrow, long ncol)
    {
        buf_.clear();
        buf_.reserve(size_t(nrow) * size_t(ncol));
        nrow_ = 0; ncol_ = ncol;
    }
    void appendRow(const T *src)
    {
        buf_.insert(buf_.end(), src, src + ncol_);
        nrow_++;
    }
    //Call once after the last appendRow: the row index cannot be built earlier
    //because insert() may reallocate the buffer and move it.
    void finishRows() { repoint(); }

    bool empty() const { return nrow_ == 0; }
    long nrow() const  { return nrow_; }
    long ncol() const  { return ncol_; }

    T       *operator[](long i)       { return row_[i]; }
    const T *operator[](long i) const { return row_[i]; }

    //For the few places that want the flat block (memcpy, gzwrite).
    T       *data()       { return buf_.empty() ? NULL : &buf_[0]; }
    const T *data() const { return buf_.empty() ? NULL : &buf_[0]; }

private:
    void repoint()
    {
        row_.resize(nrow_);
        for (long i = 0; i < nrow_; i++) row_[i] = &buf_[size_t(i) * size_t(ncol_)];
    }

    std::vector<T> buf_;
    std::vector<T *> row_;
    long nrow_;
    long ncol_;
};

#endif

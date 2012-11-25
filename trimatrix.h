#ifndef TRIMATRIX_H
#define TRIMATRIX_H

#include <QList>

template <typename T>
class TriMatrix
{
private:
    T * data;
    int matSize;
    int dataSize;

    int pos(int i, int j) const
    {
        Q_ASSERT_X(i >= 0 && i < matSize && j >= 0 && j < matSize, "TriMatrix<T>", "index out of range");
        return ( i > j ? j + i * (i + 1) / 2 : i + j * (j + 1) / 2 );
    }

public:
    TriMatrix(int size)
    {
        matSize = size;
        dataSize = size * (size + 1) / 2;
        data = new T [dataSize];
    }

    ~TriMatrix()
    {
        delete [] data;
    }

    int size() const
    {
        return matSize;
    }

    void setValue(int i, int j, const T & val) const
    {
        data[pos(i, j)] = val;
    }

    inline const T & at(int i, int j) const
    {
        return data[pos(i, j)];
    }

    void fill(const T & val) const
    {
        for (int i = 0; i < dataSize; ++i)
            data[i] = val;
    }
};

#endif // TRIMATRIX_H

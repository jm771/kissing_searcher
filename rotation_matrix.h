#include <stdint.h>
#include "types.h"
#include "vectors.h"

template <size_t Dim>
class RotationMatrix
{
    public:
    Vector<Dim> Multiply (Vector<Dim> const & vect)
    {
        Vector<Dim> ret{};

        for (size_t j = 0; j < Dim; j++) {
            for (size_t i = 0; i < Dim; i++)
            {
                ret.mValues[j] += ValueAt(i, j) * vect.mValues[i];
            }
        }

        return ret;
    }

    PointType & ValueAt(size_t i, size_t j)
    {
        return mValues[j * Dim + i];
    }

    std::array<PointType, Dim * Dim> mValues;
};

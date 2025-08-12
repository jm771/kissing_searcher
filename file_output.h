#pragma once

#include "types.h"
#include "debug_output.h"
#include <fstream>



class FileOutput
{
    public:
    FileOutput(std::string const & path) : mOutFile(path), mFirst(true), mClosed(false) {
        
    }

    template <size_t Dim>
    void WriteRow(std::vector<Vector<Dim>> const & row) {
        if (mFirst) {
            mOutFile << "[\n[";
            mFirst = false;
        } else {
            mOutFile << "],\n[";
        }
        bool first = true;
        for (auto const & vect : row) {
            if (first) {
                first=false;
            } else {
                mOutFile << ", ";
            }
            PrintVect(vect, mOutFile);

        }
    }

    void Close()
    {
        if (!mClosed) {
        mOutFile << "]\n]\n";
        mClosed = true;
        mOutFile.close();
        }
        
    }

    ~FileOutput() {
        Close();
    }


    std::ofstream mOutFile;
    bool mFirst;
    bool mClosed;
};

class NoOutput
{
    template <size_t Dim>
    void WriteRow(std::vector<Vector<Dim>> const & row) {
    }

    void Close()
    {        
    }
};
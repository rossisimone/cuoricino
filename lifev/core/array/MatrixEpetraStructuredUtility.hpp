//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
   @file MatrixEpetraStructuredUtility.hpp
   @brief The file contains utility functions to manipulate MatrixEpetraStructuredView objects

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @date 2010-11-08

   @todo createDiagBlock() and createInvDiagBlock() can be reduced to a single function to avoid copy duplication with the aid of a bit of template meta programming. In fact, a lot of routines in this file bring back to a common block and a specialized work on the single line.
 */

#ifndef _MATRIXEPETRASTRUCTUREDUTILITY_HPP_
#define _MATRIXEPETRASTRUCTUREDUTILITY_HPP_

#include <lifev/core/array/MatrixEpetraStructuredView.hpp>

namespace LifeV {

namespace MatrixEpetraStructuredUtility {

//! Copy the block specified in another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType>
void copyBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                 const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    // Processor informations
    const int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    const int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    int        srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    const int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    const int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if(( static_cast<UInt>(srcRowElement)>=srcBlock.firstRowIndex())
           && ( static_cast<UInt>(srcRowElement)<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            std::vector<int> destIndices(numSrcEntries);
            std::vector<DataType> destValues(numSrcEntries);
            int numDestEntries(0);
            int destRow(srcRowElement+rowsOffset);

            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if(( static_cast<UInt>(srcGlobalIndex)>=srcBlock.firstColumnIndex())
                   && ( static_cast<UInt>(srcGlobalIndex)<=srcBlock.lastColumnIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
            {
#ifdef NDEBUG
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
#else
                Int errorCode = destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
                ASSERT(errorCode >=0, " Error in block copy: insertion failed");
#endif
            }
            else
            {
#ifdef NDEBUG
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
#else
                Int errorCode = destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
                ASSERT(errorCode >=0, " Error in block copy: sum failed");
#endif
            }
        }
    }
}

//! Create a block full of zeros
/*!
  @param destBlock Block where the data will be stored
*/
template< typename DataType >
void createZeroBlock ( MatrixEpetraStructuredView<DataType>& destBlock )
{
    // This method will maybe be replaced
    // by the method setBlockToZero
	ASSERT( false, "The method is not yet implemented");
}

//! Create a block with an identical value on the diagonal
/*!
  @param destBlock Block where the data will be stored
  @param diagonalValue Value to be inserted in the diagonal
*/
template< typename DataType >
void createScalarBlock ( const MatrixEpetraStructuredView<DataType>& destBlock, const DataType& diagonalValue )
{
    // SQUARE TEST
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK PTR TEST
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    int destIndex(0);

    int indexBase(0);

    int firstRowIndex(destBlock.firstRowIndex()+indexBase);
    int lastRowIndex(destBlock.lastRowIndex()+indexBase);
    int firstColumnIndex(destBlock.firstColumnIndex()+indexBase);

    // Processor informations
    int  numDestElements    = destBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* destGlobalElements = destBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    int  destRowElement(0);

    for(int i(0);i<numDestElements;++i)
    {
        destRowElement = destGlobalElements[i];

        // Test if the rows are in the block
        if((destRowElement>=firstRowIndex) && (destRowElement<=lastRowIndex))
        {
            destIndex = firstColumnIndex+destRowElement-firstRowIndex;
            destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRowElement,1,&diagonalValue,&destIndex);
        }
    }
}

//! Create a block with ones on the diagonal
/*!
  @param destBlock Block where the data will be stored
*/
template< typename DataType >
void createIdentityBlock ( const MatrixEpetraStructuredView<DataType>& destBlock )
{
    createScalarBlock(destBlock,1.0);
}

//! Copy the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createDiagBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                       const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            unsigned int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            int destRow = destBlock.firstRowIndex()+diagIndex;
            int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType diagValue = 0.0;

            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is on the diagonal of the source block
                if(srcGlobalIndex-srcBlock.firstColumnIndex()==diagIndex)
                {
                    diagValue = srcValues[j];
                    j=numSrcEntries; //Exit the loop
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }
}

//! Copy the inverse of the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvDiagBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                          const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            unsigned int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            int destRow = destBlock.firstRowIndex()+diagIndex;
            int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType diagValue = 0.0;

            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is on the diagonal of the source block
                if(srcGlobalIndex-srcBlock.firstColumnIndex()==diagIndex)
                {
                    // ZERO ON DIAGONAL TEST
                	ASSERT( srcValues[j] != 0, "You cannot ask for inverse diagonal block when there are zeros on the diagonal" );

                    diagValue = 1./srcValues[j];
                    j=numSrcEntries; //Exit the loop
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }

}

//! Copy the inverse of the square root of the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvSquaredDiagBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
								 const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int  srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    unsigned int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            unsigned int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            int destRow = destBlock.firstRowIndex()+diagIndex;
            int destIndex = destBlock.firstColumnIndex()+diagIndex;
            double diagValue = 0.0;

            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is on the diagonal of the source block
                if(srcGlobalIndex-srcBlock.firstColumnIndex()==diagIndex)
                {
                    // ZERO ON DIAGONAL TEST
                	ASSERT( srcValues[j] != 0, "You cannot ask for inverse squared diagonal block when there are zeros on the diagonal" );

                    diagValue = 1/sqrt(srcValues[j]);
                    j=numSrcEntries; //Exit the loop
                }
            }
            if(srcBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }
}

//! Copy the upper part of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createUpperTriangularBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                                  const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    unsigned int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            std::vector<int> destIndices(numSrcEntries);
            std::vector<DataType> destValues(numSrcEntries);
            int numDestEntries(0);
            int destRow(srcRowElement+rowsOffset);
            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is:
                // a) in the block, and
                // b) on the upper triangular part of the block.
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()) &&
                   (srcGlobalIndex-srcBlock.firstColumnIndex()>=srcRowElement-srcBlock.firstRowIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
        }
    }
}

//! Copy the lower part of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createLowerTriangularBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                                  const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    unsigned int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            std::vector<int> destIndices(numSrcEntries);
            std::vector<DataType> destValues(numSrcEntries);
            int numDestEntries(0);
            int destRow(srcRowElement+rowsOffset);
            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is:
                // a) in the block, and
                // b) on the lower triangular part of the block.
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()) &&
                   (srcGlobalIndex-srcBlock.firstColumnIndex()<=srcRowElement-srcBlock.firstRowIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
        }
    }
}


//! Copy the lumped version of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createLumpedBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                         const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    unsigned int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            int destRow = destBlock.firstRowIndex()+diagIndex;
            int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType srcBlockRowSum = 0.0;
            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()+indexBase) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()+indexBase))
                {
                    srcBlockRowSum += abs(srcValues[j]);
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
        }
    }
}

//! Copy the inverse of the lumped version of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvLumpedBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                            const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columnss" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    unsigned int srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    DataType* srcValues;
    int* srcIndices;
    unsigned int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            int destRow = destBlock.firstRowIndex()+diagIndex;
            int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType srcBlockRowSum = 0.0;
            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()+indexBase) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()+indexBase))
                {
                    srcBlockRowSum += abs(srcValues[j]);
                }
            }

            // ZERO ON DIAGONAL TEST
        	ASSERT( srcBlockRowSum != 0, "You cannot ask for inverse lumped block when there are rows of zeros" );

            srcBlockRowSum = 1./srcBlockRowSum;
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
        }
    }
}


} // namespace MatrixEpetraStructuredUtility

} // namespace LifeV

#endif /*_MATRIXEPETRASTRUCTUREDUTILITY_HPP_ */

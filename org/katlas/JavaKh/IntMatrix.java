package org.katlas.JavaKh;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.katlas.JavaKh.algebra.rings.Int;
import org.katlas.JavaKh.interfaces.LCCC;

public class IntMatrix {
	int rows, columns;
	BigInteger matrix[][];
	IntMatrix prev, next;
	List<Integer> source, target;

	public IntMatrix(int r, int c) {
		rows = r;
		columns = c;
		matrix = new BigInteger[r][c];
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++)
				matrix[i][j] = BigInteger.ZERO;
	}

	public IntMatrix(CobMatrix<Int> cm) {
		rows = cm.target.n;
		columns = cm.source.n;
		matrix = new BigInteger[rows][columns];

		for (int i = 0; i < rows; i++) {
			LCCC<Int> rowi[] = cm.unpackRow(i);
			for (int j = 0; j < columns; j++)
				if (rowi[j] == null || rowi[j].numberOfTerms() == 0)
					matrix[i][j] = BigInteger.ZERO;
				else {
					assert rowi[j].numberOfTerms() == 1;
					matrix[i][j] = ((Int) rowi[j].firstCoefficient()).getN();
				}
		}
	}

	public boolean isDiagonal() {
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++)
				if (!matrix[i][j].equals(BigInteger.ZERO) && i != j)
					return false;
		return true;
	}

	public boolean isZero() {
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++)
				if (!matrix[i][j].equals(BigInteger.ZERO))
					return false;
		return true;
	}

	public void swapRows(int a, int b) {
		swapRows2(a, b);
		if (next != null)
			next.swapColumns2(a, b);
		if (target != null) {
			int tmp = target.get(a);
			target.set(a, target.get(b));
			target.set(b, tmp);
		}
	}

	private void swapRows2(int a, int b) {
		BigInteger tmp[] = matrix[a];
		matrix[a] = matrix[b];
		matrix[b] = tmp;
	}

	public void swapColumns(int a, int b) {
		swapColumns2(a, b);
		if (prev != null)
			prev.swapRows2(a, b);
		if (source != null) {
			int tmp = source.get(a);
			source.set(a, source.get(b));
			source.set(b, tmp);
		}
	}

	private void swapColumns2(int a, int b) {
		for (int i = 0; i < rows; i++) {
			BigInteger tmp = matrix[i][a];
			matrix[i][a] = matrix[i][b];
			matrix[i][b] = tmp;
		}
	}

	public void addRow(int a, int b, BigInteger n) {
		addRow2(a, b, n);
		/*
		 * if (next != null) --- not needed for our purposes next.addColumn2(a,
		 * b, n);
		 */
		// assert target == null || target.get(a) == target.get(b);
	}

	private void addRow2(int a, int b, BigInteger n) { // a += b * n
		for (int i = 0; i < columns; i++)
			matrix[a][i] = matrix[a][i].add(matrix[b][i].multiply(n));
	}

	public void addColumn(int a, int b, BigInteger n) {
		// assert source == null || source.get(a) == source.get(b);
		addColumn2(a, b, n);
		/*
		 * if (prev != null) --- same as above prev.addRow2(a, b, n);
		 */
	}

	private void addColumn2(int a, int b, BigInteger n) { // a += b * n
		for (int i = 0; i < rows; i++)
			matrix[i][a] = matrix[i][a].add(matrix[i][b].multiply(n));
	}

	public void multRow(int a, BigInteger n) {
		multRow2(a, n);
		/*
		 * if (next != null) next.multColumn2(a, n);
		 */
	}

	private void multRow2(int a, BigInteger n) { // a *= n
		for (int i = 0; i < columns; i++)
			matrix[a][i] = matrix[a][i].multiply(n);
	}

	public void multColumn(int a, BigInteger n) {
		multColumn2(a, n);
		/*
		 * if (prev != null) prev.multRow2(a, n);
		 */
	}

	private void multColumn2(int a, BigInteger n) { // a *= n
		for (int i = 0; i < rows; i++)
			matrix[i][a] = matrix[i][a].multiply(n);
	}

	public int rowNonZeros(int i) {
		int ret = 0;
		for (int j = 0; j < columns; j++)
			if (!matrix[i][j].equals(BigInteger.ZERO))
				ret++;
		return ret;
	}

	public int columnNonZeros(int i) {
		int ret = 0;
		for (int j = 0; j < rows; j++)
			if (!matrix[j][i].equals(BigInteger.ZERO))
				ret++;
		return ret;
	}

	public int zeroRowsToEnd() {
		int nzrows = rows;
		for (int i = 0; i < nzrows; i++)
			while (rowNonZeros(i) == 0 && i < nzrows)
				swapRows(i, --nzrows);
		return nzrows;
	}

	public int zeroColumnsToEnd() {
		int nzcols = columns;
		for (int i = 0; i < nzcols; i++)
			while (columnNonZeros(i) == 0 && i < nzcols)
				swapColumns(i, --nzcols);
		return nzcols;
	}

	public void toSmithForm(Komplex k, int param) {
		// first send all zero rows, columns to the end
		// int nzrows = zeroRowsToEnd();
		// int nzcols = zeroColumnsToEnd(); // nz means "non zero"
		// int n = Math.min(nzrows, nzcols);
		// int n = Math.min(rows, columns);
		System.out.print("\rSmith col:"+param+"           ");
		for (int row = 0, col = 0; row < rows && col < columns; row++, col++) {
            System.out.print("\rSmith col:"+param+ "\t"+row+"/"+rows+"," + col + "/"+columns+"                                         ");
			while (row < rows && rowNonZeros(row) == 0)
				row++;
			while (col < columns && columnNonZeros(col) == 0)
				col++;
			if (row >= rows || col >= columns)
				break;
			if (row > col) {
				swapRows(row, col);
				row = col;
			} else if (col > row) {
				swapColumns(row, col);
				col = row;
			}
			while (rowNonZeros(row) != 1 || columnNonZeros(col) != 1
					|| matrix[row][col].compareTo(BigInteger.ZERO) <= 0) {
				for (int j = row; j < rows; j++)
					if (matrix[j][col].compareTo(BigInteger.ZERO) < 0)
						multRow(j, BigInteger.valueOf(-1));
				while (columnNonZeros(col) != 1
						|| matrix[row][col].compareTo(BigInteger.ZERO) == 0) {
					// find the min
					BigInteger min = BigInteger.valueOf(-1);
					int idxmin = -1;
					for (int j = row; j < rows; j++)
						if ((matrix[j][col].compareTo(min) < 0 || min
								.equals(BigInteger.valueOf(-1)))
								&& matrix[j][col].compareTo(BigInteger.ZERO) > 0) {
							min = matrix[j][col];
							idxmin = j;
						}
					if (idxmin != row)
						swapRows(row, idxmin);
					// subtract the min as much as possible from the others
					for (int j = row + 1; j < rows; j++)
						if (!matrix[j][col].equals(BigInteger.ZERO))
							addRow(j, row, matrix[j][col].divide(min).negate());
				}
				// nzrows = zeroRowsToEnd();
				// nzcols = zeroColumnsToEnd();
				// n = Math.min(nzrows, nzcols);
				for (int j = col + 1; j < columns; j++)
					if (matrix[row][j].compareTo(BigInteger.ZERO) < 0)
						multColumn(j, BigInteger.valueOf(-1));
				while (rowNonZeros(row) != 1
						|| matrix[row][col].equals(BigInteger.ZERO)) {
					BigInteger min = BigInteger.valueOf(-1);
					int idxmin = -1;
					for (int j = col; j < columns; j++)
						if ((matrix[row][j].compareTo(min) < 0 || min
								.equals(BigInteger.valueOf(-1)))
								&& matrix[row][j].compareTo(BigInteger.ZERO) > 0) {
							min = matrix[row][j];
							idxmin = j;
						}
					if (idxmin != col)
						swapColumns(col, idxmin);
					for (int j = col + 1; j < columns; j++)
						if (!matrix[row][j].equals(BigInteger.ZERO))
							addColumn(j, col, matrix[row][j].divide(min)
									.negate());
				}
				// nzrows = zeroRowsToEnd();
				// nzcols = zeroColumnsToEnd();
				// n = Math.min(nzrows, nzcols);
			}
			//if ( col ==22  ) {
            if(k!=null)
                k.SaveIntMatState(param, "_" + new Int(col).toString());
//                System.out.println("Smith col:"+param+ "\t"+row+"/"+rows+"," + col + "/"+columns);
            //}

		}
		zeroRowsToEnd();
		zeroColumnsToEnd();
	}



	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				sb.append(matrix[i][j] + " ");
				sb.append("\r\n");
			}
		}
		return sb.toString();
	}

    public ArrayList<CptParam> pArray = new ArrayList<>();
	public Object synObj = null;
	public int nDone = 0;

    public void toSmithFormFast(Komplex kom, int param, int state){
        IntMatrix normMat = new IntMatrix(rows, columns);

        synObj = new Object();

        for (int row = state, col = state; row < rows && col < columns; row++, col++) {
            System.out.print("\rSmith col:"+param+ "\t"+row+"/"+rows+"," + col + "/"
                    +columns+"                                         ");
            //计算矩阵的行列 norm
//            for (int i = row; i < rows; i++) {
//                for (int j = col; j < columns; j++) {
//                    BigInteger colNorm = BigInteger.ZERO;
//                    BigInteger rowNorm = BigInteger.ZERO;
//                    for (int k = row; k < rows; k++) {
//                        colNorm = colNorm.add(matrix[k][j].multiply(matrix[k][j]));
//                    }
//                    for (int k = col; k < columns; k++) {
//                        rowNorm = rowNorm.add(matrix[i][k].multiply(matrix[i][k]));
//                    }
//                    normMat.matrix[i][j] = colNorm.multiply(rowNorm);
//                }
//            }
            //computeNormMat(normMat, row, rows, col, columns);
            //Vector<Integer> lines = new Vector<>();
            for (int i = row; i < rows; i++) {
                synchronized (synObj) {
                    pArray.add(new CptParam(i, i + 1, col, columns, row, col));
                    nDone = 0;
                }

            }

            for (int i = 0; i < kom.maxThreads; i++) {
                kom.executorService.execute(new ComputeMatNorm(this, normMat));
            }

            //ComputeMatNorm cpt = new ComputeMatNorm(this, normMat, row, rows, col, columns);
            //cpt.run();
            //选择 pivot
            synchronized (synObj) {
                try {
                    while (nDone != kom.maxThreads)
                        synObj.wait();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            BigInteger minNorm = BigInteger.ZERO;
            int pivot_i = -1;
            int pivot_j = -1;
            boolean breakFlag = false;
            //先随便选择一个非零元
			for (int i = row; i < rows && !breakFlag; i++) {
				for (int j = col; j < columns; j++) {
				    if (!matrix[i][j].equals(BigInteger.ZERO)){
				        minNorm = normMat.matrix[i][j];
				        pivot_i = i;
				        pivot_j = j;
				        breakFlag = true;
				        break;
                    }
				}
			}
			if (breakFlag == false)
			    return;
            //筛选最小的minNorm
            for (int i = row; i < rows; i++) {
                for (int j = col; j < columns; j++) {
                   if (normMat.matrix[i][j].equals(BigInteger.ZERO) || matrix[i][j].equals(BigInteger.ZERO))
                       continue;
                   if (normMat.matrix[i][j].compareTo(minNorm) < 0 ||
                           (normMat.matrix[i][j].compareTo(minNorm) == 0 &&
                                   matrix[i][j].abs().compareTo(matrix[pivot_i][pivot_j].abs()) < 0)
                           ) {
                        pivot_i = i;
                        pivot_j = j;
                        minNorm = normMat.matrix[i][j];
                   }
                }
            }
            if (pivot_i != row)
                swapRows(pivot_i, row);
            if (pivot_j != col)
                swapColumns(pivot_j, col);
            if (matrix[row][col].compareTo(BigInteger.ZERO) < 0)
                multRow(row, BigInteger.valueOf(-1));
            boolean firstTime = true;
            //迭代除法消去
            while (rowNonZeros(row) != 1 || columnNonZeros(col) !=1) {
                if (!firstTime){
                    //从所在行和列中重新选取最小 minNorm
                    //计算矩阵的行列 norm
                    synchronized (synObj) {
                        pArray.add(new CptParam(row+1, rows, col, col+1, row, col));
                        pArray.add(new CptParam(row, row+1, col+1, columns, row, col));
                        nDone = 0;
                    }

                    for (int i = 0; i < 2; i++) {
                        kom.executorService.execute(new ComputeMatNorm(this, normMat));
                    }
                    synchronized (synObj) {
                        try {
                            while (nDone != 2)
                                synObj.wait();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    }

                    //ComputeMatNorm cpt1 = new ComputeMatNorm( this, normMat, row+1, rows, col, col+1);
                    //cpt1.run();
                    //computeNormMat(normMat, row+1, rows, col, col+1);
                    //computeNormMat(normMat, row, row+1, col+1, columns);
                    //ComputeMatNorm cpt2 = new ComputeMatNorm(this, normMat, row, row+1, col+1, columns);
                    //cpt2.run();
//                    for (int i = row; i <rows ; i++) {
//                        BigInteger colNorm = BigInteger.ZERO;
//                        BigInteger rowNorm = BigInteger.ZERO;
//                        for (int k = row; k < rows; k++) {
//                            colNorm = colNorm.add(matrix[k][col].multiply(matrix[k][col]));
//                        }
//                        for (int k = col; k < columns; k++) {
//                            rowNorm = rowNorm.add(matrix[i][k].multiply(matrix[i][k]));
//                        }
//                        normMat.matrix[i][col] = colNorm.multiply(rowNorm);
//                    }
//                    for (int j = col+1; j < columns; j++) {
//                        BigInteger colNorm = BigInteger.ZERO;
//                        BigInteger rowNorm = BigInteger.ZERO;
//                        for (int k = row; k < rows; k++) {
//                            colNorm = colNorm.add(matrix[k][j].multiply(matrix[k][j]));
//                        }
//                        for (int k = col; k < columns; k++) {
//                            rowNorm = rowNorm.add(matrix[row][k].multiply(matrix[row][k]));
//                        }
//                        normMat.matrix[row][j] = colNorm.multiply(rowNorm);
//                    }
                    if (rowNonZeros(row) != 1) {
                        for (int j = col +1; j < columns; j++) {
                            if( !matrix[row][j].equals(BigInteger.ZERO)) {
                                minNorm = normMat.matrix[row][j];
                                pivot_i = row;
                                pivot_j = j;
                                break;
                            }

                        }
                    } else {
                        for (int i = row +1; i < rows; i++) {
                            if( !matrix[i][col].equals(BigInteger.ZERO)) {
                                minNorm = normMat.matrix[i][col];
                                pivot_i = i;
                                pivot_j = col;
                                break;
                            }

                        }
                    }
                    //筛选最小的minNorm
                    for (int i = row + 1; i < rows; i++) {
                        if (normMat.matrix[i][col].equals(BigInteger.ZERO) || matrix[i][col].equals(BigInteger.ZERO))
                            continue;
                        if (normMat.matrix[i][col].compareTo(minNorm) < 0 ||
                                (normMat.matrix[i][col].compareTo(minNorm) == 0 &&
                                        matrix[i][col].abs().compareTo(matrix[pivot_i][pivot_j].abs()) < 0)
                                ) {
                            pivot_i = i;
                            pivot_j = col;
                            minNorm = normMat.matrix[pivot_i][pivot_j];
                        }
                    }
                    for (int j = col + 1; j < columns; j++) {
                        if (normMat.matrix[row][j].equals(BigInteger.ZERO) || matrix[row][j].equals(BigInteger.ZERO))
                            continue;
                        if (normMat.matrix[row][j].compareTo(minNorm) < 0 ||
                                (normMat.matrix[row][j].compareTo(minNorm) == 0 &&
                                        matrix[row][j].abs().compareTo(matrix[pivot_i][pivot_j].abs()) < 0)
                                ) {
                            pivot_i = row;
                            pivot_j = j;
                            minNorm = normMat.matrix[pivot_i][pivot_j];
                        }
                    }

                    if (pivot_i != row)
                        swapRows(pivot_i, row);
                    if (pivot_j != col)
                        swapColumns(pivot_j, col);
                    if (matrix[row][col].compareTo(BigInteger.ZERO) < 0)
                        multRow(row, BigInteger.valueOf(-1));
                } else {
                    firstTime = false;
                }
                BigInteger q = matrix[col][row];
                for (int i = row + 1; i < rows; i++) {
                    if (!matrix[i][col].equals(BigInteger.ZERO)) {
                        addRow(i, row, matrix[i][col].divide(q).negate());
                        if (matrix[i][col].abs().compareTo(q.abs().divide(BigInteger.valueOf(2))) > 0) {
                            if (matrix[i][col].compareTo(BigInteger.ZERO) > 0)
                                addRow(i, row, BigInteger.valueOf(-1));
                            else
                                addRow(i, row, BigInteger.valueOf(1));
                        }
                    }
                }
                for (int j = col + 1; j < columns; j++) {
                    if (!matrix[col][j].equals(BigInteger.ZERO)) {
                        addColumn(j, col, matrix[row][j].divide(q).negate());
                        if (matrix[row][j].abs().compareTo(q.abs().divide(BigInteger.valueOf(2))) > 0) {
                            if (matrix[row][j].compareTo(BigInteger.ZERO) > 0)
                                addColumn(j, col, BigInteger.valueOf(-1));
                            else
                                addColumn(j, col, BigInteger.valueOf(1));
                        }
                    }
                }
            }

            if(kom!=null)
                kom.SaveIntMatState(param, "_" + new Int(col).toString());
        }
    }

	public static void main(String[] args) {
		IntMatrix m = new IntMatrix(2,2);
		m.matrix[0][0] = new Int(1).getN();
        m.matrix[0][1] = new Int(2).getN();
        m.matrix[1][0] = new Int(3).getN();
        m.matrix[1][1] = new Int(4).getN();

        IntMatrix n = new IntMatrix(m.rows, m.columns);
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                n.matrix[i][j] = m.matrix[i][j];
            }
        }
        m.toSmithFormFast(null, 0, 0);
        n.toSmithForm(null, 0);
        m.toString();

//        BigInteger n = new Int(-3).getN();
//        BigInteger p = new Int(-5).getN();
//        BigInteger r = p.divideAndRemainder(n)[1];
	}
}
class CptParam{
    public int istart;
    public int iend;
    public int jstart;
    public int jend;
    public int curRow;
    public int curCol;

    public CptParam(int istart, int iend, int jstart, int jend, int curRow, int curCol) {
        this.istart = istart;
        this.iend = iend;
        this.jstart = jstart;
        this.jend = jend;
        this.curRow = curRow;
        this.curCol = curCol;
    }
}
class ComputeMatNorm implements Runnable{
    private IntMatrix m;
    private IntMatrix m_norm;
    private Object synObj;


    public ComputeMatNorm(IntMatrix m, IntMatrix m_norm) {
        this.m = m;
        this.m_norm = m_norm;
        this.synObj = m.synObj;
    }


    private void computeNormMat(){
        CptParam p;
        while (true) {
            synchronized (synObj) {
                if (!m.pArray.isEmpty()) {
                    p = m.pArray.remove(0);
                } else {
                    m.nDone++;
                    synObj.notify();
                    return;
                }
            }

            for (int i = p.istart; i < p.iend; i++) {
                for (int j = p.jstart; j < p.jend; j++) {
                    BigInteger colNorm = BigInteger.ZERO;
                    BigInteger rowNorm = BigInteger.ZERO;
                    for (int k = p.curRow; k < m.rows; k++) {
                        colNorm = colNorm.add(m.matrix[k][j].multiply(m.matrix[k][j]));
                    }
                    for (int k = p.curCol; k < m.columns; k++) {
                        rowNorm = rowNorm.add(m.matrix[i][k].multiply(m.matrix[i][k]));
                    }
                    m_norm.matrix[i][j] = colNorm.multiply(rowNorm);
                }
            }
        }
    }

    @Override
    public void run() {
        computeNormMat();
    }
}
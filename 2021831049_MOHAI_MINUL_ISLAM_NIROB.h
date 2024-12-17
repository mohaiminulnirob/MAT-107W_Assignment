/*
   NAME: MOHAI MINUL ISLAM NIROB
   REG:  2021831049
*/
//MATRIX OPERATIONS (SUM,SUBTRACTION,MULTIPLICATION,DETERMINANT,COFACTOR,ADJOINT,INVERSE,SOLVING AUGMENTED MATRIX)

void printMatrix(vector<vector<double>>& matrix) {
    for (auto& row : matrix) {
        for (double element : row) {
            cout << element << "\t";
        }
        cout << "\n";
    }
}


void calculatesum(vector<vector<double>>& matrix, vector<vector<double>>&matrix1){
    int i,j;
    int m=matrix.size();
    int n=matrix[0].size();
    vector<vector<double>> sum(m, vector<double>(n, 0.0));
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        sum[i][j]=matrix[i][j]+matrix1[i][j];
    }
    cout << "\nsum matrix:\n";
    printMatrix(sum);
}


void calculatesubtraction(vector<vector<double>>& matrix, vector<vector<double>>&matrix1) {
    int i,j;
    int m=matrix.size();
    int n=matrix[0].size();
    vector<vector<double>> subtraction(m, vector<double>(n, 0.0));
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        subtraction[i][j]=matrix[i][j]-matrix1[i][j];
    }
    cout << "\nSubtraction matrix:\n";
    printMatrix(subtraction);
}


void calculatemultiplication(vector<vector<double>>& matrix, vector<vector<double>>&matrix1)
{
    int i,j,k;
    double sum=0;
    int m=matrix.size();
    int n=matrix1[0].size();
    int p=matrix[0].size();
    int q=matrix1.size();
    if(p!=q)
        cout<<"Multiplication not possible\n";
    else{
    vector<vector<double>> multiplication(m, vector<double>(n, 0.0));
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            sum=0;
            for(k=0;k<p;k++)
                sum+=matrix[i][k]*matrix1[k][j];
            multiplication[i][j]=sum;
        }
    }
    cout << "\nMultiplication matrix:\n";
    printMatrix(multiplication);
    }
}


double calculateDeterminant(vector<vector<double>>& matrix) {
    int n = matrix.size();
    if (n == 1) {
        return matrix[0][0];
    }
    double det = 0.0;
    for (int col = 0; col < n; ++col) {
        vector<vector<double>> submatrix(n - 1, vector<double>(n - 1, 0.0));

        for (int i = 1; i < n; ++i) {
            for (int j = 0, k = 0; j < n; ++j) {
                if (j != col) {
                    submatrix[i - 1][k++] = matrix[i][j];
                }
            }
        }

        det += (col % 2 == 0 ? 1 : -1) * matrix[0][col] * calculateDeterminant(submatrix);
    }
    return det;
}


void calculateCofactor(vector<vector<double>>& matrix) {
    double  det = calculateDeterminant(matrix);
    vector<vector<double>> cofactor(matrix.size(), vector<double>(matrix[0].size(), 0.0));
    if(det==0)
      cout<<"cofactor matrix doesn't exist"<<endl;
    else{
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            vector<vector<double>> submatrix(n- 1, vector<double>(n - 1, 0.0));

            for (int x = 0, m = 0; x < n; ++x) {
                if (x != i) {
                    for (int y = 0, k = 0; y < n; ++y) {
                        if (y != j) {
                            submatrix[m][k++] = matrix[x][y];
                        }
                    }
                    ++m;
                }
            }
            cofactor[i][j] =( (i + j) % 2 == 0 ? 1 : -1) * calculateDeterminant(submatrix);
        }
    }
    }
    cout << "\nCofactor matrix:\n";
    printMatrix(cofactor);
    matrix=cofactor;
}


void adjointMatrix(vector<vector<double>>&matrix) {
    double  det = calculateDeterminant(matrix);
    if(det==0)
    cout<<"transpose matrix doesn't exist"<<endl;
    else{
    calculateCofactor(matrix);
    vector<vector<double>> transpose(matrix.size(), vector<double>(matrix[0].size(), 0.0));
    int n = matrix.size();
    int m=matrix[0].size();

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            transpose[i][j] = matrix[j][i];
        }
    }
    cout << "\nAdjoint matrix:\n";
    printMatrix(transpose);
    matrix=transpose;
  }
}


void calculateInverse(vector<vector<double>>& matrix) {
    vector<vector<double>> inverse(matrix.size(), vector<double>(matrix[0].size(), 0.0));
    double det = calculateDeterminant(matrix);
    if (det == 0)
    cout << "inverse matrix doesn't exist.\n";
    else{
    adjointMatrix(matrix);
    double invDet = 1.0 / det;
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            inverse[i][j] = matrix[i][j] * invDet;
        }
    }
  }
   cout << "\nInverse matrix:\n";
   printMatrix(inverse);
   matrix=inverse;
}

void gaussianElimination(vector<vector<double>>&matrix) {
    int numRows = matrix.size();
    int numCols = matrix[0].size() - 1;

    for (int i = 0; i < numRows; i++) {
        int pivotRow = i;
        while (pivotRow < numRows && matrix[pivotRow][i] == 0) {
            pivotRow++;
        }

        if (pivotRow < numRows) {
            swap(matrix[i], matrix[pivotRow]);

            double divisor = matrix[i][i];
            for (int j = 0; j <= numCols; j++) {
                matrix[i][j] /= divisor;
            }
            for (int k = 0; k < numRows; k++) {
                if (k != i) {
                    double factor = matrix[k][i];
                    for (int j = 0; j <= numCols; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
        }
    }
    cout << "\nSolution:\n";
    for (const auto& row : matrix) {
        cout << row.back() << ' ';
    }
}

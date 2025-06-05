function A=lexicSort(A)

%LEXICSORT sorts a matrix A in lexicographic order. The code operates along
%the rows. It is based on column counterparts shared as
%http://people.sc.fsu.edu/~jburkardt/m_src/sphere_delaunay/i4col_sort_a.m
%http://people.sc.fsu.edu/~jburkardt/m_src/sphere_delaunay/i4col_compare.m
%http://people.sc.fsu.edu/~jburkardt/m_src/sphere_delaunay/i4col_swap.m
%*****************************************************************************80
%
%% I4COL_SORT_A ascending sorts an I4COL.
%
%  Discussion:
%
%    In lexicographic order, the statement "X < Y", applied to two real
%    vectors X and Y of length M, means that there is some index I, with
%    1 <= I <= M, with the property that
%
%      X(J) = Y(J) for J < I,
%    and
%      X(I) < Y(I).
%
%    In other words, the first time they differ, X is smaller.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 February 2005
%
%  Author:
%
%    John Burkardt
%   A=LEXICSORT(A)
%   * A is the matrix to be sorted
%   * A is the sorted matrix
%

assert(ismatrix(A),'The method only works for matricial arrays, while this array has %d dimensions',numDims(A));
N=size(A);
if N(1)==1;return;end
A=A';
m=N(2);n=N(1);

%FROM HERE ON WE LEAVE THE CODE ALMOST AS IT WAS
%Initialize.
i=0;indx=0;isgn=0;j=0;
%Call the external heap sorter.
while 1
    [indx,i,j]=sort_heap_external(n,indx,isgn);
    %Interchange the I and J objects.
    if indx>0;A=i4col_swap(m,n,A,i,j);
    %Compare the I and J objects.%
    elseif indx<0;isgn=i4col_compare(m,n,A,i,j);
    else break;
    end
end
A=A';

function [indx,i,j]=sort_heap_external(n,indx,isgn)

%*****************************************************************************80
%
%% SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
%
%  Discussion:
%
%    The actual list of data is not passed to the routine.  Hence this
%    routine may be used to sort integers, reals, numbers, names,
%    dates, shoe sizes, and so on.  After each call, the routine asks
%    the user to compare or interchange two items, until a special
%    return value signals that the sorting is completed.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 February 2010
%
%  Author:
%
%    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
%    MATLAB version by John Burkardt
%
%  Reference:
%
%    Albert Nijenhuis, Herbert Wilf.
%    Combinatorial Algorithms,
%    Academic Press, 1978, second edition,
%    ISBN 0-12-519260-6.
%
%  Parameters:
%
%    Input, integer N, the number of items to be sorted.
%
%    Input, integer INDX, the main communication signal.
%    The user must set INDX to 0 before the first call.
%    Thereafter, the user should set the input value of INDX
%    to the output value from the previous call.
%
%    Input, integer ISGN, results of comparison of elements I and J.
%    (Used only when the previous call returned INDX less than 0).
%    ISGN <= 0 means I is less than or equal to J;
%    0 <= ISGN means I is greater than or equal to J.
%
%    Output, integer INDX, the main communication signal.
%    If INDX is
%
%      greater than 0, the user should:
%      * interchange items I and J;
%      * call again.
%
%      less than 0, the user should:
%      * compare items I and J;
%      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
%      * call again.
%
%      equal to 0, the sorting is done.
%
%    Output, integer I, J, the indices of two items.
%    On return with INDX positive, elements I and J should be interchanged.
%    On return with INDX negative, elements I and J should be compared, and
%    the result reported in ISGN on the next call.
%

    persistent i_save;
    persistent j_save;
    persistent k;
    persistent k1;
    persistent n1;

    if isempty(i_save);i_save=-1;end
    if isempty(j_save);j_save=-1;end
        
    if indx==0%INDX = 0: This is the first call.  
        k = floor ( n / 2 );
        k1 = k;
        n1 = n;
    elseif indx<0%INDX < 0: The user is returning the results of a comparison.
        if indx==-2
            if isgn<0;i_save=i_save+1;end
            j_save=k1;
            k1=i_save;
            indx=-1;
            i=i_save;
            j=j_save;
            return;
        end
        if 0<isgn
            indx=2;
            i=i_save;
            j=j_save;
            return;
        end
        if k<=1
            if n1==1
                i_save=0;
                j_save=0;
                indx=0;
            else
                i_save=n1;
                n1=n1-1;
                j_save=1;
                indx=1;
            end
            i=i_save;
            j=j_save;
            return;
        end
        k=k-1;
        k1=k;
    elseif ( indx == 1 )%0 < INDX, the user was asked to make an interchange.
        k1 = k;
    end
    while 1
        i_save=2*k1;
        if i_save==n1
            j_save=k1;
            k1=i_save;
            indx=-1;
            i=i_save;
            j=j_save;
            return;
        elseif i_save<=n1
            j_save=i_save+1;
            indx=-2;
            i=i_save;
            j=j_save;
            return;
        end
        if k<=1;break;end
        k=k-1;
        k1=k;
    end
    if n1==1
        i_save=0;
        j_save=0;
        indx=0;
        i=i_save;
        j=j_save;
    else
        i_save=n1;
        n1=n1-1;
        j_save=1;
        indx=1;
        i=i_save;
        j=j_save;
    end
end

function a=i4col_swap(m,n,a,i,j)

%*****************************************************************************80
%
%% I4COL_SWAP swaps columns I and J of a integer array of column data.
%
%  Example:
%
%    Input:
%
%      M = 3, N = 4, I = 2, J = 4
%
%      A = (
%        1  2  3  4
%        5  6  7  8
%        9 10 11 12 )
%
%    Output:
%
%      A = (
%        1  4  3  2
%        5  8  7  6
%        9 12 11 10 )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns in the array.
%
%    Input, integer A(M,N), an array of N columns of length M.
%
%    Input, integer I, J, the columns to be swapped.
%
%    Output, integer A(M,N), the array, with columns I and J swapped.
%

    if i<1 || n<i || j<1 || n<j
        fprintf(1,'\n');
        fprintf(1,'I4COL_SWAP - Fatal error!\n');
        fprintf(1,'  I or J is out of bounds.\n');
        fprintf(1,'  I =    %d\n',i);
        fprintf(1,'  J =    %d\n',j);
        fprintf(1,'  N =    %d\n',n);
        error('I4COL_SWAP - Fatal error!');
    end

    if i==j;return;end

    col(1:m) = a(1:m,i)';
    a(1:m,i) = a(1:m,j);
    a(1:m,j) = col(1:m)';
end

function isgn=i4col_compare(m,n,a,i,j)

%*****************************************************************************80
%
%% I4COL_COMPARE compares columns I and J of a integer array.
%
%  Example:
%
%    Input:
%
%      M = 3, N = 4, I = 2, J = 4
%
%      A = (
%        1  2  3  4
%        5  6  7  8
%        9 10 11 12 )
%
%    Output:
%
%      ISGN = -1
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 June 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns.
%
%    Input, integer A(M,N), an array of N columns of vectors of length M.
%
%    Input, integer I, J, the columns to be compared.
%    I and J must be between 1 and N.
%
%    Output, integer ISGN, the results of the comparison:
%    -1, column I < column J,
%     0, column I = column J,
%    +1, column J < column I.
%

%
%  Check.
%
    if i<1
        fprintf(1,'\n');
        fprintf(1,'I4COL_COMPARE - Fatal error!\n');
        fprintf(1,'  Column index I = %d < 1.\n',i);
        error('I4COL_COMPARE - Fatal error!');
    end
    if n<i
        fprintf(1,'\n');
        fprintf(1,'I4COL_COMPARE - Fatal error!\n');
        fprintf(1,'  N = %d < column index I = %d.\n',n,i);
        error('I4COL_COMPARE - Fatal error!');
    end

    if j<1
        fprintf(1,'\n');
        fprintf(1,'I4COL_COMPARE - Fatal error!\n');
        fprintf(1,'  Column index J = %d < 1.\n',j);
        error('I4COL_COMPARE - Fatal error!');
    end

    if n<j
        fprintf(1,'\n');
        fprintf(1,'I4COL_COMPARE - Fatal error!\n');
        fprintf(1,'  N = %d < column index J = %d.\n',n,j);
        error('I4COL_COMPARE - Fatal error!');
    end
    isgn=0;
    if i==j;return;end
    k=1;

    while k<=m
        if a(k,i)<a(k,j);isgn=-1;return
        elseif a(k,j)<a(k,i);isgn=+1;return
        end
        k=k+1;
    end
end

end
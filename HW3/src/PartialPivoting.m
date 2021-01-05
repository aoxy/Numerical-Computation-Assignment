clear, clc

A = [1,1,1;1,1,2;1,2,2];
b = [1;2;1];
showAb(A,b);
% Implement Gauss elimination with partial pivoting
[sx,sy] = size(A);
P = eye(sx);
% Create a long string of dashes and prepare to print to the command
% window.
dashString = repmat('-', 1, 62);
for j = 1:sy - 1
    [~,row] = max(abs(A(j:sx,j)));
    row = row + j - 1;
    U = A(row,:);
    if row ~= j
        fprintf('Row %d and %d to be interchanged.\n',j,row);
        fprintf('%s\n', dashString);
        pause;
        A(row,:) = A(j,:);
        A(j,:) = U;
        V = P(row,:);
        P(row,:) = P(j,:);
        P(j,:) = V;
        t = b(row);
        b(row) = b(j);
        b(j) = t;
        fprintf('Row interchange:\n');
        showAb(A,b);
        fprintf('%s\n', dashString);
    end
    pause;
    fprintf('Elimination:\n');
    for i = j + 1:sx
        if A(i,j) == 0
            continue
        end
        multiplier = - A(i,j)/A(j,j);
        A(i,j) = -multiplier;
        A(i,j + 1:sy) = A(i,j + 1:sy) + multiplier*A(j,j + 1:sy);
        b(i) = b(i) + multiplier*b(j);
    end
    showAb(A,b);
end
pause;
L = tril(A,-1) + diag(ones(sx,1));
U = triu(A);
fprintf('\n%s\n', dashString);
fprintf('A = P''LU\n')
showLU(L,U,P);

function showLU(L,U,varargin)
% Display the matrix L and U got in LU decomposition
if ~isempty(varargin)
    P = varargin{1};
else
    P = -1;
end
% Get the size
[sx,sy] = size(L);
for i = 1:sx
    if i == round((1 + sx)/2)
        Lstring = 'L = ';
        Ustring = ', U = ';
        if P ~= -1
            Pstring = ', P = ';
        end
    else
        Lstring = '    ';
        Ustring = '      ';
        if P ~= -1
            Pstring = '      ';
        end
    end
    for j = 1:sy
        nextstring = sprintf('%c',rats(L(i,j)));
        t = find(nextstring ~= ' ');
        nextstring = nextstring(t(1):t(end));
        if length(nextstring) < 5
            fillstring1 = repmat(' ',1,floor((5 - length(nextstring))/2));
            fillstring2 = repmat(' ',1,ceil((5 - length(nextstring))/2));
            nextstring = [fillstring1,nextstring,fillstring2];
        end
        Lstring = [Lstring, nextstring];
        nextstring2 = sprintf('%c',rats(U(i,j)));
        tt = find(nextstring2 ~= ' ');
        nextstring2 = nextstring2(tt(1):tt(end));
        if length(nextstring2) < 5
            fillstring1 = repmat(' ',1,floor((5 - length(nextstring2))/2));
            fillstring2 = repmat(' ',1,ceil((5 - length(nextstring2))/2));
            nextstring2 = [fillstring1,nextstring2,fillstring2];
        end
        Ustring = [Ustring,nextstring2];
        if P ~= -1
            nextstring3 = sprintf('%c',rats(P(i,j)));
            tt = find(nextstring3 ~= ' ');
            nextstring3 = nextstring3(tt(1):tt(end));
            if length(nextstring3) < 5
                fillstring1 = repmat(' ',1,floor((5 - length(nextstring3))/2));
                fillstring2 = repmat(' ',1,ceil((5 - length(nextstring3))/2));
                nextstring3 = [fillstring1,nextstring3,fillstring2];
            end
            Pstring = [Pstring,nextstring3];
        end
    end
    if P ~= -1
        laststring = [Lstring,'  ',Ustring,'  ',Pstring];
    else
        laststring = [Lstring,'  ',Ustring];
    end
    % Print to the command window:
    fprintf(laststring);
    fprintf('\n');
end
fprintf('\n');
end

function showAb(A,b)
% Display the matrix A and the vector b used in Gauss elimination

% Get the size
[sx,sy] = size(A);
for i = 1:sx
    if i == round((1 + sx)/2)
        showstring = 'A = ';
        bstring = ', b = ';
    else
        showstring = '    ';
        bstring = '      ';
    end
    for j = 1:sy
        nextstring = sprintf('%c',rats(A(i,j)));
        t = find(nextstring ~= ' ');
        nextstring = nextstring(t(1):t(end));
        if length(nextstring) < 5
            fillstring1 = repmat(' ',1,floor((5 - length(nextstring))/2));
            fillstring2 = repmat(' ',1,ceil((5 - length(nextstring))/2));
            nextstring = [fillstring1,nextstring,fillstring2];
        end
        showstring = [showstring, nextstring];
    end
    nextstring2 = sprintf('%c',rats(b(i)));
    tt = find(nextstring2 ~= ' ');
    nextstring2 = nextstring2(tt(1):tt(end));
    if length(nextstring2) < 5
        fillstring1 = repmat(' ',1,floor((5 - length(nextstring2))/2));
        fillstring2 = repmat(' ',1,ceil((5 - length(nextstring2))/2));
        nextstring2 = [fillstring1,nextstring2,fillstring2];
    end
    bstring = [bstring,nextstring2];
    laststring = [showstring,'  ',bstring];
    % Print to the command window:
    fprintf(laststring);
    fprintf('\n');
end
fprintf('\n');
end
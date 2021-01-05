clear, clc

A = [2,1,1;2,2,-1;4,-1,6];
b = [9;9;16];
showAb(A,b);
% Implement Gauss elimination without pivoting or storing the multipliers.
[sx,sy] = size(A);
for j = 1:sy - 1
    pause
    fprintf('GE on column %d: \n',j);
    for i = j + 1:sx
        if A(i,j) == 0
            continue
        end
        multiplier = - A(i,j)/A(j,j);
        A(i,j) = 0;
        A(i,j + 1:sy) = A(i,j + 1:sy) + multiplier*A(j,j + 1:sy);
        b(i) = b(i) + multiplier*b(j);
    end
    showAb(A,b);
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
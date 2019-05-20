function p = fft_split(x,M)

N = 4^M;
L=length(x);

if(N<L)
    error('N should be  greater than L ');
end
%to ensure that the subsequent splittings occur without any problem
xn = x;
xn=[xn zeros(1,N-L)];
xe=xn(2:2:end);
xo1=xn(1:4:end);
xo3=xn(3:4:end);

%splitting the input matrix to sub parts even and odd(Cooley Tukey
%algorithm)
for k=0:N-1
      for r=0:N/2-1
          Wne=exp(-1j*pi*k*r*4/N);
          Xe(k+1,r+1) = Wne;
      end
end
%subsequent splittings based on N value

for k=0:N-1
      for r=0:N/4-1
          Wno=exp(-1j*pi*k*r*8/N);
          Xo(k+1,r+1) = Wno;
      end
end

for k=0:N-1
    Wn1=exp(-1j*2*pi*k/N);
    Wn3=exp(-1j*6*pi*k/N);
end
    
%computing the exponential factor (using roots of unity concept)
p = Xe*xe' + Wn1*(Xo*xo1')+Wn3*(Xo*xo3');
%final fft computation
end

subroutine newton(x)
implicit none
real x,x1,f,f1
integer::i=0
f(x)=4*cos(x)-exp(x)
f1(x)=-4*sin(x)-exp(x)
x1=x-f(x)/f1(x)
do while(abs(x-x1)>1.0e-6.and.i<=1000)
    x=x1
    i=i+1
    x1=x-f(x)/f1(x)
enddo
if(i<=20)then
print*,'某点附近根x=',x1
else
    print*,'附近无根'
endif
end subroutine
program ntroot
implicit none
real x0
print*,'请输入某点x'
read*,x0
call newton(x0)
end program ntroot

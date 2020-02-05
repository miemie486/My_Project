plot_point1.cpp is the program to calculate (x-1)^9 and expansion of (x-1)^9， in single, double, quad. 
plot.py is the program to plot those point which was calculated by plot_point1.cpp.
plot_point_Horner_Method.cpp is the problem use Horner's Method to calculate the expansion of (x-1)^9 in single.
plot_Horner.py is the program to plot those point which was calculated by plot_point_Horner_Method.cpp. And this program also compare it with those point in datafloat1.txt and datafloat2.txt.

Problem 2:
(3)
I find that the error of Horner's Method in single are much larger than the expansion of (x-1)^9 in single. If I subtract two number in computer, there will have a subtractive cancelation. And Hoener's Method will multiply each subtractive cancelation with x^n (n depends on the place of subtractive cancelation.)

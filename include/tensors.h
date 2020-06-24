// Useful Tensors

#define dab(a,b) (a==b)?1:0

double eabc(int a,int b,int c){
	double reteab;
	if((a==XUP&&b==YUP&&c==ZUP)||(a==YUP&&b==ZUP&&c==XUP)||(a==ZUP&&b==XUP&&c==YUP))
		reteab=1.0;
	else if((a==XUP&&b==ZUP&&c==YUP)||(a==YUP&&b==XUP&&c==ZUP)||(a==ZUP&&b==YUP&&c==XUP))
		reteab=-1.0;
	else
		reteab=0;
	return reteab;
}
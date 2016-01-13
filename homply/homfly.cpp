 

/* Program written by H.R.Morton and H.B.Short,
   University of Liverpool, May 1985.
 Amended most recently, March 1997, minor modification Nov 1997. 
 Translated to C++ 2002 , Feb 2003 and modified */
 #include <iostream>
 #include <string>
 #include <iomanip>
 #include<algorithm>
using namespace std;
const int maxl = 160; // longest braid string to be input
const int Vassmax=maxl;  // Truncation degree for Vassiliev invariants
const int globaln = 9;  // largest number of strings used 
const int globalnfact = 362880; // =globaln!

struct polynomial
{
				long coefft[Vassmax+1];
				int maxdeg;
				int mindeg;
};                             
                          


int n,len, r;
int writhe;
int vassiliev_type=Vassmax;
bool mistake, calculated[globaln+1];
char ch;

int braid[maxl];
string braidstr;
int mult[globalnfact][globaln];
int factor[globaln+1];


polynomial c[globalnfact];
polynomial zeropoly; 
polynomial u[globaln+1];

int w[globaln+maxl][globaln];
int results[globaln+maxl][globaln];
int conw[globaln+maxl];




int main()
{
   int i,j,k,s;
   int bin(int,int);
   
   void multiply(int,int);
   void wrapup();
   void init();
   void braidinput();
   void reinit();
   
   void screenout();
   void multtab();
   polynomial minus(polynomial);
   polynomial plus(polynomial,polynomial);
   polynomial shift(polynomial);
   
   n=globaln;
   
   init();
   
   for( ; ; )
   {
       mistake=false;
       braidinput(); // mistake =true if braid input is wrong
       if (mistake==false)
       { 
          reinit();
          cout  << "The polynomial of the closed braid is being calculated" << endl;
          for( i = 1; i <= len; i++) 
          { 
              r = braid[i];
              multiply(1,factor[n]);                                     
          } // end for( i = 1; i <= len; i++) 
   
          wrapup(); // This closes off the braid strings.
   
          for( j = 0; j <= Vassmax; j++) 
          {          
              for( s = 0; s <= n-1; s++) 
              {
                  k = n-s-1;                                               
                  if (u[s].coefft[j] != 0)  
                  {
                      for( i = 0; i <= k; i++)
                      {
                          w[j+s][i] += u[s].coefft[j]* bin(k,i);
                      } // end for( i = 0; i <= k; i++)
                  } // end if (u[s].coefft[j] != 0)
              } // end for( s = 0; s <= n-1; s++)
          } //end for( j = 0; j <= Vassmax; j++)
          for( i = 0; i <= Vassmax; i++) 
          {
              conw[i] = 0;  
              for( j = 0; j <= n-1; j++) 
              {    
                   results[i][j] = w[i][j]; 
                   conw[i]+= w[i][j];
              }// end for( j = 0; j <= n-1; j++)
          } // end for( i = 0; i <= Vassmax; i++)
      
          screenout();
       } // end if (mistake==false)
       cout << "Do you want to calculate another example? Enter y/n." << endl;
       cin >>ch;
       if (ch=='n')  break;                 
   
   } // end for( ; ; )


}// end of main program

void init()                                                                 
{
   int i,j ;
   for (j=1; j<= globaln; j++)   calculated[j]= false; 
   // calculated keeps track of how far the multiplication table mult has been found
 
   zeropoly.mindeg=-1000;zeropoly.maxdeg=-1000;
   for (i=0; i<= Vassmax; i++)  zeropoly.coefft[i]=0; // def of zeropoly

   factor[0]=1;                                                               
   for (j= 1; j<= n; j++)  factor[j]= j*factor[j-1] ;  
   // definition of factorial           

} //end {init}  

void braidinput()
{  
   int i,j;
   char ans, ch;
   int sig=1;
   cout <<endl; 
   cout << "This program calculates the Homfly polynomial of a closed braid" << endl;
   cout << "on at most " <<globaln <<" strings."<< endl;
   cout << endl;                                                     

   cout << "Give braid to be calculated, for example 213-241-3, " << endl;
   cout << "using 1,2,.. for the standard braid generators" << endl;
   cout << "and -1,-2,... for their inverses." << endl;   
   cout << endl;
  
   cin >> braidstr;                 
  
   cout << "The number of strings to be used in this braid is currently " << n << endl;
   cout << "If you want to continue with this number type y" << endl;
   cout << "Otherwise give the number of strings to be used in the braid" << endl;
         
   cin >> ans; 
   if (ans!='y') 
   {  
      n=ans-'0';
      if ((n>globaln) || (n<2)) 
      {
         cout  << "The number of strings should be >1 and at most " << globaln << endl;
         mistake=true;
      } // end     if ((n>globaln) || (n<2))    
   }// end  if (ans!='y') 
  
   len =braidstr.size();          
   writhe = 0; 
   i = 0; 
   for( j =0 ; j <= len-1; j++) 
   { 
       ch = braidstr[j]; 
       if (ch == '-')  sig = -1;
       else 
       { 
          i++;
          braid[i] = sig*(ch - '0');
          writhe += sig; 
          sig = 1;                                       
          if (ch <='0'||ch - '0' > n-1)  
          {
              cout  << "wrong element "<< ch << " in braid at position " <<  
              j<<endl; mistake=true;  
          } // end if (ch <='0'||ch - '0' > n-1) 
       } // end if (ch == '-') 
    } // end  for( j = 0; j <= len-1; j++) 
    len = i; 
    // len is now the actual braid length, not the string length which has - signs

    cout << endl;
    if (mistake==false) cout << "Braid length is " << len << endl;
    cout << endl;
    if (len > maxl-15)  cout  << "DANGER..braid may be too long!!" << endl;         
} // end braidinput

void reinit()
{   
   int i,j;
   void multtab();   
        
   for (j=0; j<= n; j++)  u[j]=zeropoly; // initialise array u of polys

   for(j=0;j<=globaln+maxl-1;j++)
   {
     for (i=0;i<=globaln-1;i++)
       {
        w[j][i]=0;
       }// end for (i=0;i<=globaln-1;i++)
    }// end for(j=0;j<=globaln+maxl-1;j++)
    // initialise working array w of coefficients
                                          
   for (i=1; i<= factor[n]-1;i++)   c[i]=zeropoly;                                              
   c[0].mindeg=0;
   c[0].coefft[0]=1;
   c[0].maxdeg=0;
   for (j=1; j<= Vassmax; j++)  c[0].coefft[j]=0;
   // initialise main working array c of n! polynomials

 
   cout<< "Please wait while multiplication table is calculated" << endl;
   multtab();     // calculation of multiplication table mult                                                               

   cout << endl;
   cout << endl;
}// end reinit

int bin(int v, int v1)
{
   int j,denom,num;
   int bin_result;
   denom = 1;
   num = 1;
   for( j = 1; j <= v; j++) num = num*j;
   for( j = 1; j <= v1; j++) denom = -denom*j;
   for( j = 1; j <= v-v1; j++) denom = denom*j;
   bin_result=  num / denom;
   return bin_result;
}//end bin (binomial coefficient)



void multtab()                                                             
{
   int cr,cr1,m,g,f,ff,h,j,k;                                                                                                                  
   if (calculated[n]==false)
   {                        
      k=1;                                             
      for (g= 0 ;g<= factor[n]-1; g++)
      {                                                                        
         mult[g][1]=g+k;                                                           
         k=-k;                                                                    
     
       }//end    for (g= 0 ;g<= factor[n]-1; g++)
                                                                        
      for (r=2; r<= n-1; r++)
      {                                                                        
         f=factor[r]; ff =factor[r-1];                                           
         for (g= 0; g<= factor[n]-1; g++)
         {                                                                    
            m= g / f;                                                       
            cr= m % (r+1);                                                      
            m= g / ff;                                                      
            cr1= m % r;                                                         
            j= cr1 - cr;                                                          
            if (j>= 0)
            {                                                                 
               h= g + j*ff*(r-1) + f;                                             
               mult[g][r]= h;                                                      
               mult[h][r]= g;                                                      
             }//end   if(j>=0) 
          }// end   for (g= 0; g<= factor[n]-1; g++)                                                                  
       }// end      for (r=2; r<= n-1;r++)  
    for (j=1;j<=n;j++) calculated[j]=true;//avoids recalculation of mult later
    }//end  if (calculated[n]==false)                                                          
} //end {multtab}       

void multiply(int a,int b)                                              
{
   int sign,g,h;                                                          
   polynomial  x,y;  
   polynomial minus(polynomial);
   polynomial plus(polynomial,polynomial);
   polynomial shift(polynomial);
                               
   if (r>0)  sign= 1; else sign= -1;                                         
   r=abs(r);   
                                                                  
   for (g=a-1;g<=b-1; ++g)
   {                                                                      
      h= mult[g][r];                                
      if (h > g)                                                              
      {                                                                     
         if (c[g].mindeg>-1 || c[h].mindeg>-1) 
         {                                                        
             swap(c[g],c[h]);
             if (sign>0)                                                          
             {   
                c[h]=plus(c[h],shift(c[g]));
              } else                                                               
             {             
                c[g]=plus(c[g],minus(shift(c[h])));
             }//end  if (sign>0) 
         }//end  if (c[g].mindeg>-1 || c[h].mindeg>-1)                                                                 
      }// end if {h > g}                                                                     
   } //end  for (g=a-1;g<=b-1; ++g)                                                                        
}// end {multiply}  

void wrapup()
{                                                               
   int i,j,k,s,a,b,g,rr,f,ff;                                        
   polynomial plus(polynomial,polynomial)  ;                                                                      
   for (s=n-1;s>=1;s--)                                                     
   {                                                                        
      b=0; f=factor[s];                                                       
      for (i=0; i<= n-s-1; i++)                                                      
      {                                                                    
         b=b+factor[s+i]; a=b+1-f;                                            
         for (rr=1; rr<= s-1;rr++)                                                    
         {                                                                 
            k=(s-rr+1)*f;r=rr;                                                
            multiply(a+k,b+k);                                                  
            for (g=a+k-1; g<=  b+k-1; g++) // corrected for indexing from 0 in c                                              
            {
               if (c[g].mindeg>=0) 
               c[g-f]=plus(c[g],c[g-f]);
            }// end for (g=a+k-1; g<=  b+k-1; g++)
         }// end    for (rr=1; rr<= s-1;rr++)                                                              
         ff=f-factor[s+i];                                                     
         for (g=a-1;g<= b-1; g++)  // corrected for indexing from 0 in c
         {
            if (i>0)                                      
            {
               if (c[g].mindeg>=0) 
               {
                  c[g+ff]=plus(c[g],c[g+ff]);
               }// end if c[g].mindeg>=0
            }//end if (i>0)
         }// end for (g=a-1;g<= b-1; g++) 
      }// end    for (i=0; i<= n-s-1; i++)                                                              
   }// end   for (s=n-1;s>=1;s--)                                                                        
    f=0;                                                                        
    for (j=0; j<= n-1; j++)                                                           
    {                                                                       
       f=f+factor[j];                                                           
       u[j]=c[f-1];   // corrected for indexing from 0 in c                                                         
    }// end for (j=0; j<= n-1; j++)                                                                         
}// end {wrapup}                

polynomial shift(polynomial a) //{multiplies a polynomial in z by the variable z}
{
   int i;
   polynomial x;
   if (a.mindeg<0)  return zeropoly; else
   {
      x=zeropoly;
      x.mindeg=a.mindeg+1;
      x.maxdeg=min(a.maxdeg+1,vassiliev_type);
   /* allows truncation of the invariant at a given degree in z,
    specified during calculation */
   
      for (i=x.mindeg; i<= x.maxdeg; i++) 
      x.coefft[i]=a.coefft[i-1];//end for (i=x.mindeg; i<= x.maxdeg; i++)
      return x;
   } //end if (a.mindeg<0)
}// end {shift}  

polynomial plus(polynomial a,polynomial b)                                    
{
   polynomial c; int i;
   if (a.mindeg<0)  return b; else
   {
      if (b.mindeg<0)  return a; else
      {
          c = zeropoly;
          c.mindeg=min(a.mindeg,b.mindeg);
          c.maxdeg=max(a.maxdeg,b.maxdeg);
          for (i=c.mindeg; i<= c.maxdeg; i++) 
          c.coefft[i]=a.coefft[i]+b.coefft[i]; // end for (i=c.mindeg; i<= c.maxdeg; i++) 
          return c;
      } // end if (b.mindeg<0)
   } // end if (a.mindeg<0) 
} // end {plus}

polynomial minus(polynomial a)                                    
{
   polynomial c; int i;
   if (a.mindeg<0)  return zeropoly; else
   {
      c = zeropoly;
      c.mindeg=a.mindeg;
      c.maxdeg=a.maxdeg;
      for (i=c.mindeg; i<= c.maxdeg; i++) 
      c.coefft[i]=-a.coefft[i]; // end for (i=c.mindeg; i<= c.maxdeg; i++) 
      return c;
   } // end if (a.mindeg<0) 
} // end {minus}

                                                           

void screenout()                                                 
{
   int j,k,u,s;
   bool nonzero;                                                            
   cout << endl;                                                           
   cout << "braid: " << braidstr << endl; cout  << endl;         
   cout << "writhe: " << writhe<< endl; cout  << endl; 
   cout << "Array of coefficients of the Homfly polynomial," << endl;
   cout << "with columns indexed by degree in v and rows by degree in z:" << endl;
   cout << endl;
  
   for( u=0; u <= n-1; u++) cout  <<setw(8)<< writhe-n+1+2*u;                     
   cout  << endl;
   for( u=0; u <= n-1; u++) cout  <<setw(8)<<"________";                     
   cout  << setw(3)<<"___"<< endl;            
   for( k=0; k <= Vassmax; k++)                                                       
   {
      j=0;s=1;
      for( u=0; u <= n-1;  u++)                                     
      { 
         if (results[k][u]!=0)                                                    
         {
            j=1; cout  <<setw(8*s)<<results[k][u];
            s=1;                            
         }
         else s=s+1; // end  if (results[k][u]!=0)                                                  
      } //end for( u=0; u <= n-1;  u++)                                                
      if (j==1)  cout <<setw(8*(s-1)+3) << "|" <<setw(3)<< k-n+1 << endl;                             
   }// end for( k=0; k <= Vassmax; k++)

   for( u=0; u <= n-1; u++) cout  <<setw(8)<<"________";                     
   cout  << setw(3)<<"___"<< endl; // to rule off the coefficient array
 
   cout  << endl;                                                      
   cout  << "Coefficients of Conway polynomial:" << endl; 
   cout << endl;  
   nonzero=false;                                   
   for( k=0; k <= Vassmax; k++)                                                       
   if (conw[k]!=0)
   {
       cout <<setw(8) << conw[k] << " (" << setw(3)<< k-n+1<< ")" << endl;
       nonzero=true;
   } // end if (conw[k]!=0)
   if (nonzero== false) cout << "Conway = 0" << endl;
   cout<<endl;                                                                
}  // end {screenout}                                                                          







#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
    unsigned int blue;
    unsigned int green;
    unsigned int red;
}pixel;

typedef struct
{
    int x, y, nr;
    double cor;
}detectii;

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
   FILE *fin, *fout;
   unsigned int dim_img, latime_img, inaltime_img;
   unsigned char pRGB[3], header[54], aux;

   printf("nume_fisier_sursa = %s \n",nume_fisier_sursa);

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fout = fopen(nume_fisier_destinatie, "wb+");

   fseek(fin, 2, SEEK_SET);
   fread(&dim_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

   //copiaza octet cu octet imaginea initiala in cea noua
	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	//calculam padding-ul pentru o linie
	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;

    else
        padding = 0;

    printf("padding = %d \n",padding);

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			//citesc culorile pixelului
			fread(pRGB, 3, 1, fout);
			//fac conversia in pixel gri
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
}

unsigned int* XORSHIFT32(unsigned int seed, int n)
{
    int k;
    unsigned long int *r=(unsigned long int *)calloc(n, sizeof(unsigned long int));
    r[0]=seed;
    for(k=1;k<n;k++)
    {
        seed=seed^seed<<13;
        seed=seed^seed>>17;
        seed=seed^seed<<5;
        r[k]=seed;
    }
    return r;
}

pixel** matrice(int W, int H)
{
    pixel **I=(char**)calloc(H, sizeof(pixel*));
    for(int i=0;i<H;i++)
        I[i]=(pixel*)calloc(W, sizeof(pixel));
    return I;
}

pixel** citire(FILE *f, int W, int H)
{
    pixel **I=matrice(W, H);
    fseek(f, 54, SEEK_SET);
    int i, j, k, n;
    unsigned char *aux=(unsigned char*)calloc(3, sizeof(unsigned char));
    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
        {
            fread(aux, 3, 1, f);
            I[i][j].red=aux[2];
            I[i][j].green=aux[1];
            I[i][j].blue=aux[0];
        }
    return I;
}

pixel* liniarizare(int W, int H, pixel **I)
{
    pixel *L=(pixel*)calloc(W*H, sizeof(pixel));
    int k=0;
    for(int i=H-1;i>=0;i--)
        for(int j=W-1;j>=0;j--)
        {
            L[k].blue=I[i][j].blue;
            L[k].green=I[i][j].green;
            L[k].red=I[i][j].red;
            k++;
        }
    return L;
}

int* Durstenfeld(int W, int H, unsigned long int *R)
{
    int r, k, aux;

   int *p=(int *)calloc(W*H, sizeof(int));
    for(k=0;k<W*H;k++)
        p[k]=k;

    for(k=W*H-1;k>=1;k--)
    {
        r=R[W*H-k]%(k+1);
        aux=p[r];
        p[r]=p[k];
        p[k]=aux;
    }
    return p;
}

pixel* permutare(pixel *L, int H, int W, unsigned int *R, int *p)
{
    int i;
    pixel *P=(pixel*)calloc(W*H, sizeof(pixel));
    for(i=0;i<W*H;i++)
    {
        P[p[i]]=L[i];
        P[p[i]].green=L[i].green;
        P[p[i]].red=L[i].red;
    }
    return P;
}

pixel* inversa(int W, int H, pixel *C, int *p)
{
    int i, j, aux, *pd=(int*)calloc(W*H, sizeof(int));
    pixel *Pd=(pixel*)calloc(W*H, sizeof(pixel));
    for(i=0;i<W*H;i++)
        pd[p[i]]=i;
    for(i=0;i<W*H;i++)
    {
        Pd[pd[i]].blue=C[i].blue;
        Pd[pd[i]].green=C[i].green;
        Pd[pd[i]].red=C[i].red;
    }
    return Pd;
}

pixel* criptare(pixel *P, unsigned int *R, unsigned int SV, int W, int H)
{
    pixel *C=(pixel*)calloc(W*H, sizeof(pixel));
    C[0].blue=(SV)^(P[0].blue)^(R[W*H]);
    C[0].green=(SV>>8)^(P[0].green)^(R[W*H]>>8);
    C[0].red=(SV>>16)^(P[0].red)^(R[W*H]>>16);
    for(int k=1;k<W*H;k++)
    {
        C[k].blue=(C[k-1].blue)^(P[k].blue)^(R[W*H+k]);
        C[k].green=(C[k-1].green)^(P[k].green)^(R[W*H+k]>>8);
        C[k].red=(C[k-1].red)^(P[k].red)^(R[W*H+k]>>16);
    }
    return C;
}

pixel* decriptare(pixel *C, unsigned int *R, unsigned int SV, int W, int H)
{
    pixel *D=(pixel*)calloc(W*H, sizeof(pixel));
    D[0].blue=(SV)^(C[0].blue)^(R[W*H]);
    D[0].green=(SV>>8)^(C[0].green)^(R[W*H]>>8);
    D[0].red=(SV>>16)^(C[0].red)^(R[W*H]>>16);
    for(int k=1;k<W*H;k++)
    {
        D[k].blue=(C[k-1].blue)^(C[k].blue)^(R[W*H+k]);
        D[k].green=(C[k-1].green)^(C[k].green)^(R[W*H+k]>>8);
        D[k].red=(C[k-1].red)^(C[k].red)^(R[W*H+k]>>16);
    }
    return D;
}

void creared(FILE *fin, FILE *fout, pixel *C, int H, int W)
{
    unsigned char c;
    int i;
    fseek(fin,0,SEEK_SET);
    while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);
	fseek(fout, 54, SEEK_SET);
    for(i=W*H-1;i>=0;i--)
    {
        fwrite(&C[i].blue, 1, 1, fout);
        fwrite(&C[i].green, 1, 1, fout);
        fwrite(&C[i].red, 1, 1, fout);
        fflush(fout);
    }
    fclose(fout);
}

void creare(FILE *fin, FILE *fout, pixel **C, int H, int W)
{
    unsigned char c;
    int i, j;
    fseek(fin,0,SEEK_SET);
    while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);
	fseek(fout, 54, SEEK_SET);
    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
        {
            fwrite(&C[i][j].blue, 1, 1, fout);
            fwrite(&C[i][j].green, 1, 1, fout);
            fwrite(&C[i][j].red, 1, 1, fout);
            fflush(fout);
        }
    fclose(fout);
}

pixel* frecventa(pixel *L, int W, int H)
{
    pixel *fr=(pixel*)calloc(256, sizeof(pixel));
    for(int i=0;i<256;i++)
        for(int j=0;j<W*H;j++)
        {
            if(L[j].red==i)
                fr[i].red++;
            if(L[j].green==i)
                fr[i].green++;
            if(L[j].blue==i)
                fr[i].blue++;
        }
    return fr;
}

void test_chi(int m, int n, pixel *I)
{
    double ff=(m*n)/256, chi_r=0, chi_g=0, chi_b=0;
    pixel *fr=frecventa(I, m, n);
    for(int i=0;i<256;i++)
    {
        chi_r+=((fr[i].red-ff)*(fr[i].red-ff))/ff;
        chi_g+=((fr[i].green-ff)*(fr[i].green-ff))/ff;
        chi_b+=((fr[i].blue-ff)*(fr[i].blue-ff))/ff;
    }
    printf("R: %.2lf\nG: %.2lf\nB: %.2lf\n", chi_r, chi_g, chi_b);
}

double media_intensitatiilor(int i, int j, int n, int m, pixel **I)
{
    double s=0;
    double a=n*m;
    int k, l;
    for(k=i;k<n;k++)
        for(l=j;l<m;l++)
            s+=I[k][l].red;
    s/=a;
    return s;
}

double deviatia_standard(int i, int j, int n, int m, pixel **I)
{
    double a=(n*m)-1, med;
    med=media_intensitatiilor(i, j, n, m, I);
    a=(1.0)/a;
    int k, l;
    double s=0;
    for(k=i;k<n;k++)
        for(l=j;l<m;l++)
            s+=((I[k][l].red-med)*(I[k][l].red-med));
    s*=a;
    s=sqrt(s);
    return s;
}

pixel** colorare(int i, int j, pixel **C, int c)
{
    int k, l, r, g, b;
    if(c==0)
    {
        r=255;
        g=0;
        b=0;
    }
    if(c==1)
    {
        r=255;
        g=255;
        b=0;
    }
    if(c==2)
    {
        r=0;
        g=255;
        b=0;
    }
    if(c==3)
    {
        r=0;
        g=255;
        b=255;
    }
    if(c==4)
    {
        r=255;
        g=0;
        b=255;
    }
    if(c==5)
    {
        r=0;
        g=0;
        b=255;
    }
    if(c==6)
    {
        r=192;
        g=192;
        b=192;
    }
    if(c==7)
    {
        r=255;
        g=140;
        b=0;
    }
    if(c==8)
    {
        r=128;
        g=0;
        b=128;
    }
    if(c==9)
    {
        r=128;
        g=0;
        b=0;
    }
    for(k=j;k<j+11;k++)
    {
        C[i][k].red=r;
        C[i][k].green=g;
        C[i][k].blue=b;
    }
    for(k=i;k<i+15;k++)
    {
        C[k][j+10].red=r;
        C[k][j+10].green=g;
        C[k][j+10].blue=b;
    }
    for(k=j;k<j+11;k++)
    {
        C[i+14][k].red=r;
        C[i+14][k].green=g;
        C[i+14][k].blue=b;
    }
    for(k=i;k<i+15;k++)
    {
        C[k][j].red=r;
        C[k][j].green=g;
        C[k][j].blue=b;
    }
    return C;
}

detectii* glisare(pixel **I, int H, int W, int h, int w, pixel **S, int c, detectii *D, int *k, double ps)
{
    int i, j, js, is;
    double sum;
    double corr, a=h*w-1;
    a=(1.0)/a;
    double med_s, med_f, dev_s, dev_f;
    dev_s=deviatia_standard(0, 0, h, w, S);
    med_s=media_intensitatiilor(0, 0, h, w, S);
    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
        {
            sum=0;
            dev_f=deviatia_standard(i, j, i+h, j+w, I);
            med_f=media_intensitatiilor(i, j, i+h, j+w, I);
            double b=(1.0)/(dev_s*dev_f);
            for(is=0;is<h;is++)
                for(js=0;js<w;js++)
                    sum+=(b*(I[i+is][j+js].red-med_f)*(S[is][js].red-med_s));
            sum+=(b*(I[i+is][j+js].red-med_f)*(-med_s));
            corr=a*sum;
            if(corr>=ps && corr<=1)
            {
                D[(*k)].cor=corr;
                D[(*k)].x=j;
                D[(*k)].y=i;
                D[(*k)].nr=c;
                (*k)++;
            }
        }
    return D;
}

int cmp(const void *a, const void *b)
{
    detectii va=*(detectii*)a, vb=*(detectii*)b;
    if(va.cor<vb.cor)
        return 1;
    if(va.cor>vb.cor)
        return -1;
    return 0;
}

void suprapunere(detectii **D, int *k)
{
    double sp;
    int i, j, m, n, l;
    for(i=0;i<(*k)-1;i++)
    {
        double aria_i=165;
        for(j=i+1;j<(*k);j++)
        {
            double aria_j=165;
            double aria_int, lungime, latime;
            lungime=(*D)[i].y+14-(*D)[j].y;
            latime=(*D)[i].x+10-(*D)[j].x;
            aria_int=lungime*latime;
            sp=aria_i+aria_j-aria_int;
            aria_int/=sp;
            if(aria_int>0.2)
            {
                for(l=j;l<(*k)-1;l++)
                {
                    (*D)[l].cor=(*D)[l+1].cor;
                    (*D)[l].x=(*D)[l+1].x;
                    (*D)[l].y=(*D)[l+1].y;
                    (*D)[l].nr=(*D)[l+1].nr;
                }
            (*k)--;
            }
        }
    }
}

int main()
{
    ///prima parte
    FILE *nr=fopen("secret_key.txt", "r");

    unsigned int R0, SV;
    fscanf(nr, "%lu %lu", &R0, &SV);

    int W, H;
    unsigned int *R, *Rd;

    char *img=(char*)calloc(20, sizeof(char)), img_criptata1[]="criptare.bmp", img_decriptata1[]="decriptare.bmp";
    pixel **I, *L, *P, *C, *Pd, *D, **PC, *LC;

    printf("Prima imagine: ");
    scanf("%s", img);

	FILE *fin, *fout;
	fin=fopen(img, "rb");
    fseek(fin, 18, SEEK_SET);
    fread(&W, sizeof(unsigned int), 1, fin);
    fread(&H, sizeof(unsigned int), 1, fin);

    R=(unsigned int *)calloc(32, sizeof(unsigned int));
    R=XORSHIFT32(R0, 2*W*H-1);

    int *p=Durstenfeld(W, H, R);

    I=citire(fin, W, H);

    L=liniarizare(W, H, I);
    printf("Chi-squared test on RGB channels for peppers.bmp:\n");
    test_chi(W, H, L);

    P=permutare(L, H, W, R, p);

    C=criptare(P, R, SV, W, H);

    fout=fopen(img_criptata1, "wb+");
    creared(fin, fout, C, H, W);

    ///decriptare
    FILE *fd, *gd;
    fd=fopen(img_criptata1, "rb+");

    fseek(fd, 18, SEEK_SET);
    fread(&W, sizeof(unsigned int), 1, fd);
    fread(&H, sizeof(unsigned int), 1, fd);

    Rd=(unsigned int *)calloc(32, sizeof(unsigned int));
    Rd=XORSHIFT32(R0, 2*W*H-1);

    PC=citire(fd, W, H);
    LC=liniarizare(W, H, PC);
    printf("Chi-squared test on RGB channels for criptare.bmp:\n");
    test_chi(W, H, LC);

    p=Durstenfeld(W, H, R);

    D=decriptare(C, Rd, R0, W, H);
    P=permutare(D, H, W, Rd, p);
    Pd=inversa(W, H, D, p);

    gd=fopen(img_decriptata1, "wb+");
    creared(fd, gd, Pd, H, W);

    ///partea a 2a
    int H2, W2;
    char *img_sursa=(char*)calloc(20,sizeof(char)), img_dest[]="test_grayscale.bmp", img_criptata[]="criptare2.bmp", img_decriptata[]="decriptare2.bmp";
    char *cif0=(char*)calloc(20,sizeof(char)), *cif1=(char*)calloc(20,sizeof(char)), *cif2=(char*)calloc(20,sizeof(char)), *cif3=(char*)calloc(20,sizeof(char)), *cif4=(char*)calloc(20,sizeof(char)), *cif5=(char*)calloc(20,sizeof(char)), *cif6=(char*)calloc(20,sizeof(char)), *cif7=(char*)calloc(20,sizeof(char)), *cif8=(char*)calloc(20,sizeof(char)), *cif9=(char*)calloc(20,sizeof(char));
    char cif0d[]="cifra0d.bmp", cif1d[]="cifra1d.bmp", cif2d[]="cifra2d.bmp", cif3d[]="cifra3d.bmp", cif4d[]="cifra4d.bmp", cif5d[]="cifra5d.bmp", cif6d[]="cifra6d.bmp", cif7d[]="cifra7d.bmp", cif8d[]="cifra8d.bmp", cif9d[]="cifra9d.bmp";

    printf("Imaginea sursa: ");
    scanf("%s", img_sursa);

    unsigned int *R2;
    R2=(unsigned int *)calloc(32, sizeof(unsigned int));
    R2=XORSHIFT32(R0, 2*W*H-1);

    FILE *f, *g, *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9, *fdd;
    pixel **I2, **S0, **S1, **S2, **S3, **S4, **S5, **S6, **S7, **S8, **S9, **C2, *L2, **GS;
	f=fopen(img_sursa, "rb");
    fseek(f, 18, SEEK_SET);
    fread(&W2, sizeof(unsigned int), 1, f);
    fread(&H2, sizeof(unsigned int), 1, f);
    I2=citire(f, W, H);


    ///criptare


    I2=citire(f, W2, H2);


    L2=liniarizare(W2, H2, I2);
    printf("Chi-squared test on RGB channels for test.bmp:\n");
    test_chi(W2, H2, L2);

    pixel *P2;

    int *p2=Durstenfeld(W2, H2, R2);

    P2=permutare(L2, H2, W2, R2, p2);

    C2=criptare(P2, R2, SV, W2, H2);

    g=fopen(img_criptata, "wb+");
    creared(f, g, C2, H2, W2);

    ///decriptare
    FILE *gdd;
    fdd=fopen(img_criptata, "rb+");

    fseek(fdd, 18, SEEK_SET);
    fread(&W2, sizeof(unsigned int), 1, fdd);
    fread(&H2, sizeof(unsigned int), 1, fdd);

    unsigned int *Rd2;
    Rd2=(unsigned int *)calloc(32, sizeof(unsigned int));
    Rd2=XORSHIFT32(R0, 2*W2*H2-1);

    pixel **PC2, *LC2, *D2, *Pd2;
    PC2=citire(fdd, W2, H2);
    LC2=liniarizare(W2, H2, PC2);
    printf("Chi-squared test on RGB channels for criptare.bmp:\n");
    test_chi(W2, H2, LC2);

    D2=decriptare(C2, Rd2, SV, W2, H2);
    Pd2=inversa(W2, H2, D2, p2);

    gdd=fopen(img_decriptata, "wb+");
    creared(fdd, gdd, Pd2, H2, W2);


    ///template matching
    printf("Sablonul 0: ");
    scanf("%s", cif0);

    printf("Sablonul 1: ");
    scanf("%s", cif1);

    printf("Sablonul 2: ");
    scanf("%s", cif2);

    printf("Sablonul 3: ");
    scanf("%s", cif3);

    printf("Sablonul 4: ");
    scanf("%s", cif4);

    printf("Sablonul 5: ");
    scanf("%s", cif5);

    printf("Sablonul 6: ");
    scanf("%s", cif6);

    printf("Sablonul 7: ");
    scanf("%s", cif7);

    printf("Sablonul 8: ");
    scanf("%s", cif8);

    printf("Sablonul 9: ");
    scanf("%s", cif9);

    grayscale_image(img_decriptata, img_dest);
    grayscale_image(cif0, cif0d);
    grayscale_image(cif1, cif1d);
    grayscale_image(cif2, cif2d);
    grayscale_image(cif3, cif3d);
    grayscale_image(cif4, cif4d);
    grayscale_image(cif5, cif5d);
    grayscale_image(cif6, cif6d);
    grayscale_image(cif7, cif7d);
    grayscale_image(cif8, cif8d);
    grayscale_image(cif9, cif9d);

    f0=fopen(cif0d, "rb+");
    fseek(f0, 18, SEEK_SET);
    int w, h;
    fread(&w, sizeof(unsigned int), 1, f0);
    fread(&h, sizeof(unsigned int), 1, f0);

    S0=citire(f0, w, h);

    f1=fopen(cif1d, "rb+");
    S1=citire(f1, w, h);

    f2=fopen(cif2d, "rb+");
    S2=citire(f2, w, h);

    f3=fopen(cif3d, "rb+");
    S3=citire(f3, w, h);

    f4=fopen(cif4d, "rb+");
    S4=citire(f4, w, h);

    f5=fopen(cif5d, "rb+");
    S5=citire(f5, w, h);

    f6=fopen(cif6d, "rb+");
    S6=citire(f6, w, h);

    f7=fopen(cif7d, "rb+");
    S7=citire(f7, w, h);

    f8=fopen(cif8d, "rb+");
    S8=citire(f8, w, h);

    f9=fopen(cif9d, "rb+");
    S9=citire(f9, w, h);

    FILE *fdest=fopen(img_dest, "rb+"), *gdest=fopen("glisare.bmp", "wb+");

    pixel **Idest;

    Idest=citire(fdest, W2, H2);

    GS=citire(fdest, W2, H2);


    int k=0;
    detectii *Dd;


    Dd=(detectii*)calloc(999999, sizeof(detectii));

    Dd=glisare(Idest, H2-h, W2-w, h, w, S0, 0, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S1, 1, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S2, 2, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S3, 3, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S4, 4, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S5, 5, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S6, 6, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S7, 7, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S8, 8, Dd, &k, 0.50);

    Dd=glisare(Idest, H2-h, W2-w, h, w, S9, 9, Dd, &k, 0.50);

    qsort(Dd, k, sizeof(detectii), cmp);

    suprapunere(&Dd, &k);
    printf("%d", k);
    Dd=realloc(Dd, k*sizeof(detectii));

    for(int i=0;i<k;i++)
        GS=colorare(Dd[i].y, Dd[i].x, GS , Dd[i].nr);

    FILE *e=fopen(img_decriptata, "rb+");
    creare(e, gdest, GS, H2, W2);

    return 0;
}

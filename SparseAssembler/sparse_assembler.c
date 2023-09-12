#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "redblack_t.h"

char linija[1024];

int snimout;
int percout;
double mtx_min, mtx_max, mtx_epsilon;
int mtx_count = 0;

int bSymetric;

/* Spurse Structures */
int nCounter;
rb_callback_func_t rb_callback_func;
RBNode_t *root;
RBNode_t *nil;
RBNode_t nilObj;

int * packRow;
int * packCol;
double * packVal;

char *SparseFileName = "sparse.bin";
FILE *SparseStream;

void CountNode_callback(RBNode_t * pNode) {nCounter++;}
void BrisiNode_callback(RBNode_t * pNode) {free(pNode);}
void SetZeroVal_callback(RBNode_t * pNode) {pNode->info = 0;}

int RBCount(RBNode_t * pNode)
{
	nCounter = 0;
	rb_callback_func = CountNode_callback;
	inorder(pNode);
	return nCounter;
}

void SetAllZero()
{
	rb_callback_func = SetZeroVal_callback;
	inorder(root);
}

void RBInit()
{
	root = nil = &(nilObj);
	nil->key = -1;
	nil->color = black;
}

void clean(RBNode_t * p)
{
	rb_callback_func = BrisiNode_callback;
	postorder(p);
}

void insert(Tkey k, Tinfo inf)
{
	RBNode_t *p = (RBNode_t*)malloc(sizeof(RBNode_t));
	p->key=k;
	p->info=inf;
	rb_insert(p);
}

void insert2(RBNode_t * q)
{
	RBNode_t *p, *r;
	r = nil;
	p = root;
	while (p!=nil) {
		r = p;
		if (q->key<p->key) p = p->leftson; else p = p->rightson;
	}
	q->father = r;
	q->leftson = nil;
	q->rightson = nil;
	if (r==nil) root = q;
	else if (q->key<r->key) r->leftson = q; else r->rightson = q;
}

RBNode_t * search(RBNode_t *p, Tkey k)
{
	while ( p!=nil )
	{
		if (k==p->key) break;
		if (k<p->key) p = p->leftson; else p = p->rightson;
	}
	return p;
}

void inorder(RBNode_t *p)
{
	if (p!=nil) {
		inorder(p->leftson);
		rb_callback_func(p);
		inorder(p->rightson);
	}
}

void preorder(RBNode_t *p)
{
	if (p!=nil) {
		rb_callback_func(p);
		preorder(p->leftson);
		preorder(p->rightson);
	}
}

void postorder(RBNode_t *p)
{
	if (p!=nil) {
		postorder(p->leftson);
		postorder(p->rightson);
		rb_callback_func(p);
	}
}

RBNode_t * minimum(RBNode_t *p)
{
	while (p->leftson!=nil) p = p->leftson;
	return p;
}

RBNode_t * maximum(RBNode_t *p)
{
	while (p->rightson!=nil) p = p->rightson;
	return p;
}

RBNode_t * successor(RBNode_t *p)
{
	RBNode_t *q;
	if (p->rightson!=nil) return minimum(p->rightson);
	else {
		q = p->father;
		while ( (q!=nil) && (p==q->rightson) ) {
			p = q;
			q = q->father;
		}
		return q;
	}
}

RBNode_t * predecessor(RBNode_t *p)
{
	RBNode_t *q;
	if (p->leftson!=nil) return maximum(p->leftson);
	else {
		q = p->father;
		while ( (q!=nil) && (p==q->leftson) ) {
			p = q;
			q = q->father;
		}
		return q;
	}
}


void leftrotate(RBNode_t *x)
{
	RBNode_t *y;
	y = x->rightson;
	x->rightson = y->leftson;
	if (y->leftson!=nil) y->leftson->father = x;
	y->father = x->father;
	if (x->father==nil) root = y;
	else if (x==x->father->leftson) x->father->leftson = y;  
	else x->father->rightson = y;
	y->leftson = x;
	x->father = y;
}

void rightrotate(RBNode_t *x)
{
	RBNode_t *y;
	y = x->leftson;
	x->leftson = y->rightson;
	if (y->rightson!=nil) y->rightson->father = x;
	y->father = x->father;
	if (x->father==nil) root = y;
	else if (x==x->father->leftson) x->father->leftson = y;  
	else x->father->rightson = y;
	y->rightson = x;
	x->father = y;
}

void rb_insert(RBNode_t *x)
{
	RBNode_t *y;
	insert2(x);
	x->color = red;
	while ( (x!=(root)) && ( x->father->color==red) ) 
		if (x->father==x->father->father->leftson) {
			y = x->father->father->rightson;
			if (y->color==red) {
				x->father->color = black;        /*  Fall 1  */
				y->color = black;                /*          */
				x->father->father->color = red;  /*          */
				x =  x->father->father;                   /*          */
			}
			else {
				if (x==x->father->rightson) {
					x = x->father;                          /*  Fall 2  */
					leftrotate(x);                          /*          */
				}
				x->father->color = black;        /*  Fall 3  */ 
				x->father->father->color = red;  /*          */
				rightrotate(x->father->father);           /*          */
			}
		}
		else {
			y = x->father->father->leftson;
			if (y->color==red) {
				x->father->color = black;        /*  Fall 1  */
				y->color = black;                /*          */
				x->father->father->color = red;  /*          */
				x =  x->father->father;                   /*          */
			}
			else {
				if (x==x->father->leftson) {
					x = x->father;                          /*  Fall 2  */
					rightrotate(x);                         /*          */
				}
				x->father->color = black;        /*  Fall 3  */ 
				x->father->father->color = red;  /*          */
				leftrotate(x->father->father);            /*          */
			}
		}
		root->color = black;
}

void rb_delete(RBNode_t *z)
/* entfernt *z im Baum */
{
	RBNode_t *x, *y;
	if ( (z->leftson==nil) || (z->rightson==nil) ) y = z;
	else y = successor(z);
	if (y->leftson!=nil) x = y->leftson; else x = y->rightson;
	x->father = y->father;
	if (y->father==nil) root = x;
	else if (y==y->father->leftson) y->father->leftson = x;
	else                       y->father->rightson = x;
	if (y!=z) {
		z->key = y->key;
		z->info = y->info;
	}
	if (y->color==black) rb_delete_fixup(x);
	free(y);
}

void rb_delete_fixup(RBNode_t *x)
{
	RBNode_t *w;
	while ( (x!=root) && (x->color==black) )
		if (x==x->father->leftson) {
			w = x->father->rightson;
			if (w->color==red) {
				w->color = black;                /*  Fall 1  */
				x->father->color = red;          /*          */
				leftrotate(x->father);                    /*          */
				w = x->father->rightson;                  /*          */
			}
			if (    (w->leftson->color==black) 
				&& (w->rightson->color==black) ) {
					w->color = red;                  /*  Fall 2  */
					x = x->father;                            /*          */
			}
			else {
				if (w->rightson->color==black) {
					w->leftson->color = black;     /*  Fall 3  */
					w->color = red;                /*          */
					rightrotate(w);                         /*          */
					w = x->father->rightson;                /*          */
				}
				w->color = x->father->color;              /*  Fall 4  */
				x->father->color = black;        /*          */
				w->rightson->color = black;      /*          */
				leftrotate(x->father);                    /*          */
				x = root;
			}
		}
		else {
			w = x->father->leftson;
			if (w->color==red) {
				w->color = black;                /*  Fall 1  */
				x->father->color = red;          /*          */
				rightrotate(x->father);                   /*          */
				w = x->father->leftson;                   /*          */
			}
			if (    (w->leftson->color==black) 
				&& (w->rightson->color==black) ) {
					w->color = red;                  /*  Fall 2  */
					x = x->father;                            /*          */
			}
			else {
				if (w->leftson->color==black) {
					w->rightson->color = black;    /*  Fall 3  */
					w->color = red;                /*          */
					leftrotate(w);                          /*          */
					w = x->father->leftson;                 /*          */
				}
				w->color = x->father->color;              /*  Fall 4  */
				x->father->color = black;        /*          */
				w->leftson->color = black;       /*          */
				rightrotate(x->father);                   /*          */
				x = root;
			}
		}
		x->color = black;
}  

void SnimiElem_callback(RBNode_t * pNod)
{
	packRow[snimout]  = *(((int*)&pNod->key) + 1);
	packCol[snimout] = *((int*)&pNod->key);
	packVal[snimout] = pNod->info;
	snimout++;
}

void AddVal(int row, int col, Tinfo val)
{
	RBNode_t * pNode;
	Tkey index;
	if (val == 0) return;
	*(((int*)&index) + 1) = row;
	*(int*)&index = col;
	pNode = search(root, index);
	if (pNode->key == index)
	{
		pNode->info += val;
		if (pNode->info == 0) rb_delete(pNode);
	}
	else insert(index, val);
}

void SetVal(int row, int col, Tinfo val)
{
	RBNode_t * pNode;
	Tkey index;
	if (val == 0) return;
	*(((int*)&index) + 1) = row;
	*(int*)&index = col;
	pNode = search(root, index);
	if (pNode->key == index)
		pNode->info = val;
	else
		insert(index, val);
}

// For GFORTRAN: void __cdecl sparseassembler_kill()
void __cdecl SPARSEASSEMBLER_KILL()
{
    printf("Killing Sparse Assembler...\n");
	clean(root);
}

//void __cdecl sparseassembler_init(int *symetric)
void __cdecl SPARSEASSEMBLER_INIT(int *symetric)
{
    printf("Initializing Sparse Assembler...\n");
	RBInit();
	bSymetric = *symetric;
}

//void __cdecl sparseassembler_addelemmatrix(int *n, int *indices, double *vals)
void __cdecl SPARSEASSEMBLER_ADDELEMMATRIX(int *n, int *indices, double *vals)
{
    int i,j, nn = *n;

	for(i=0;i<nn;i++)
	{
		for(j=(bSymetric ? i : 0);j<nn;j++)
		{
            if ( indices[i]!=0 && indices[j]!= 0 )
            {
//                if ((indices[i] <= indices[j]) ||  (!bSymetric))
//    		        AddVal(indices[i], indices[j], vals[i*nn+j]);
//	            else
			        AddVal(indices[j], indices[i], vals[i*nn+j]);
            }
		}
	}
}

//void __cdecl sparseassembler_getsparse(int *nz, int *rows, int *cols, double *vals) 
void __cdecl SPARSEASSEMBLER_GETSPARSE(int *nz, int *rows, int *cols, double *vals)
{
	packRow = rows;
	packCol = cols;
	packVal = vals;
	snimout=0;
	rb_callback_func = SnimiElem_callback;
	inorder(root);
    *nz = snimout;
    printf("Nonzero count: %d\n", *nz);
}

void __cdecl SPARSEASSEMBLER_GETNONZERO(int *nz) {
	*nz = RBCount(root);
}

void __cdecl SPARSEASSEMBLER_SETVAL(int *row, int *col, double *val)
{
	SetVal(*row, *col, *val);
}

void __cdecl SPARSEASSEMBLER_ZERO()
{
	SetAllZero();
}

// Dodao Milos, snima sparse matricu u fajl zbog ustede memorije
void SnimiElem_u_Fajl_callback(RBNode_t * pNod)
{
	int row, col;
    double val;

    row  = *(((int*)&pNod->key) + 1);
	col = *((int*)&pNod->key);
	val = pNod->info;

    fwrite(&row, sizeof(row), 1, SparseStream);
    fwrite(&col, sizeof(col), 1, SparseStream);
    fwrite(&val, sizeof(val), 1, SparseStream);
}

// Dodao Milos, snima sparse matricu u fajl zbog ustede memorije
void __cdecl SPARSEASSEMBLER_SAVESPARSEFILE()
{
	rb_callback_func = SnimiElem_u_Fajl_callback;
    if( (SparseStream  = fopen( SparseFileName, "wb" )) == NULL )
    {
        printf("Cannot open %s file for writing!\n", SparseFileName);
        return;
    }
	inorder(root);
    fclose(SparseStream);
}

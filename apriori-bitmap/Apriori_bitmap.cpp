/* author: hsinhoyeh yhh92u@gmail.com
 * goal  : this cpp implement the bitmap version of apriori algorithm
 * */

#include "stdafx.h"

#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <time.h>

using namespace std;

#define default_item -1

//max number of items in each transaction
#define n_max_items_per_tran 50

//the number of bits in each int structure
#define int_size 32

//from 0 to 999
#define sz_of_item 1000

//#define minsup 0.004
double minsup = 0.01;

typedef unsigned int uint;
typedef struct _itemset{
//change the structure of items into array
	int* items;
	int* bitmap;
} Itemset;


list< int* > trans;
list< list<Itemset>* > Ln_itemsets;

//store the L1 frequence
int * L1;
int ** trans_bitmap = NULL;
int size_of_trans = -1;
int freq_for_minsup=-1;
int bitmap_length = -1;

static int pos2bitmap[] = {
	0x00000001, 
	0x00000002,
	0x00000004,
	0x00000008,
	0x00000010,
	0x00000020,
	0x00000040,
	0x00000080,
	0x00000100,
	0x00000200,
	0x00000400,
	0x00000800,
	0x00001000,
	0x00002000,
	0x00004000,
	0x00008000,
	0x00010000,
	0x00020000,
	0x00040000,
	0x00080000,
	0x00100000,
	0x00200000,
	0x00400000,
	0x00800000,
	0x01000000,
	0x02000000,
	0x04000000,
	0x08000000,
	0x10000000,
	0x20000000,
	0x40000000,
	0x80000000,
};
//release the itemset
void Itemset_free(Itemset& a){
	if(a.items != NULL)
		delete[] a.items;
	if(a.bitmap != NULL)
		free(a.bitmap);
}
//initialize the itemset
void Itemset_init(Itemset& a, int level){
	a.items = new int[level];
	for(int i=0; i<level; i++){
		a.items[i] = -1;
	}
}
//push an item into the itemset
void Itemset_push(Itemset& a, int pos, int item){
	a.items[pos] = item;
}
//push b into a
void Itemset_pushArr(Itemset& a, const Itemset& b, int level){
	for(int lcnt=0; lcnt < level ; lcnt++){
		a.items[lcnt] = b.items[lcnt];
	}
}

//scan the transaction data into mem
//and build the L1
void loadData(const char* filename){

	L1 = (int*)malloc(sizeof(int) * sz_of_item);
	for(int inicnt=0; inicnt < sz_of_item; inicnt++){
		L1[inicnt] = 0;
	}

	int* tran = NULL;
	int c;
	int cur_item = default_item;
	FILE* pFile = fopen(filename, "rb");
	if (pFile==NULL) perror ("Error opening file");
	else{
		do {
		  c = fgetc(pFile);
		  
		  if(c == '\n'){
			//new transaction meet

			//flush the buffer
			  if(cur_item != default_item){

				  if(tran == NULL){
					tran = (int*)malloc(sizeof(int)* (n_max_items_per_tran+1));
					for(int ini_item=1; ini_item<= n_max_items_per_tran; ini_item++ ){
						tran[ini_item] = default_item;
					}
					tran[0] = 1;
				  }	


				  L1[cur_item] ++;
				  tran[ tran[0] ] = cur_item;
				  cur_item = default_item;
				  tran[0] ++;
			  }


			  trans.push_back(tran);
			  tran = NULL;
		  }
		  else if(c == '\r'){
			//ignore
		  }
		  else if(c == ' '){
			//new item delima, push into current transaction
			  if(tran == NULL){
				tran = (int*)malloc(sizeof(int)* (n_max_items_per_tran+1));
				for(int ini_item=1; ini_item<= n_max_items_per_tran; ini_item++ ){
					tran[ini_item] = default_item;
				}
				tran[0] = 1;
			  }
			//push into L1
			L1[cur_item] ++;


			//push into transaction
			tran[ tran[0] ] = cur_item;
			cur_item = default_item;
			tran[0] ++;
			//the total number is less than tran[0]

			
		  }
		  else{
			//item content
			  if(cur_item == default_item){
				cur_item = c - 0x30;
			  }
			  else{
				cur_item = cur_item * 10 + (int)(c - 0x30);
			  }
		  }

		} while (c != EOF);
		
		size_of_trans = trans.size();
		freq_for_minsup = (int)(size_of_trans * (double)minsup);

		fclose (pFile);
	}
}
//release data
void freeData(){
	list< int* >::const_iterator citer = trans.begin();
	for(; citer != trans.end(); ++citer){
		free( (*citer) );
	}
	trans.clear();
}
//dump the content of transaction into file
void dumpTrans(const char* filename){
	FILE* pFile = fopen(filename, "wb");
	if (pFile==NULL) perror ("Error opening file");
	else{
		list< int* >::const_iterator citer = trans.begin();
		for(; citer != trans.end(); ++citer){
			int* tran = (*citer);
			for(int itemcnt=1; itemcnt < tran[0]; itemcnt++){
				fprintf(pFile, "%d ", tran[itemcnt]);
			}
			fprintf(pFile, "\n");
		}
		fclose(pFile);
	}
}
//bitmap init
int* bitmap_init(int bitmap_length){
	int* ret = (int*)malloc(sizeof(int) * bitmap_length);
	for(int iint=0; iint < bitmap_length; iint++){
		ret[iint] = 0;
	}
	return ret;
}
//push a transaction id into bitmap
void bitmap_push(int* bitmap, int tranid){
	if(bitmap == NULL){
		fprintf(stderr, "push bitmap is null\n");
		return ;
	}
	int index = tranid / int_size;
	int pos = tranid % int_size;
//	fprintf(stderr, "insert tranid=%d, index=%d, pos=%d\n",tranid, index, pos);
	bitmap[index] = bitmap[index] | pos2bitmap[pos];
}
//erase a transaction id from bitmap
void bitmap_pop(int* bitmap, int tranid){
	if(bitmap == NULL){
		fprintf(stderr, "push bitmap is null\n");
		return ;
	}
	int index = tranid / int_size;
	int pos = tranid % int_size;
//	fprintf(stderr, "pop tranid=%d, index=%d, pos=%d\n",tranid, index, pos);
	bitmap[index] = bitmap[index] & ( 0xffffffff ^ pos2bitmap[pos]);
}
//join two bitmap into another one
int* bitmap_intersect(int* b1, int* b2, int bitmap_length){
	int * ret_b = bitmap_init(bitmap_length);
	int freq=0;
	for(int icnt=0; icnt< bitmap_length;icnt++){
		if(b1[icnt] == 0)
			continue;
		else if(b2[icnt] == 0)
			continue;
		else{
			ret_b[icnt] = b1[icnt]&b2[icnt];	
			uint bmap = (uint)ret_b[icnt];
			while(bmap != 0){
				if((bmap & 0x00000001)){
					freq++;
				}
				bmap = bmap >> 1;
			}
		}
	}
	if(freq < freq_for_minsup){
		//in-frequent
		free(ret_b);
		return NULL;
	}
	else
		return ret_b;
}
//convert the horizental structure into vertical structure
//and generate L1
void convertData2Bitmap(){
	//convert the row data into bitmap format and allocate l1 itemset
	bitmap_length = (int)ceil((double)size_of_trans / (double)int_size);
	fprintf(stderr, "bitmap length = %d\n", bitmap_length);
	trans_bitmap = (int**)malloc(sizeof(int*) * sz_of_item);
	for(int bcnt=0; bcnt < sz_of_item; bcnt++){
		trans_bitmap[bcnt] = NULL;
	}
	
	

	int item;
	int tranid=0;
	//scan the whole data once to built the bitmap
	list< int* >::const_iterator citer = trans.begin();
	for(; citer != trans.end(); ++citer){
		int* tran = (*citer);
		for(int tcnt=1; tcnt < tran[0]; tcnt++){
			item = tran[tcnt];
			if(L1[item] < freq_for_minsup)
				continue;
				//frequent < minsup, ignore
			//built bitmap
			if(trans_bitmap[item] == NULL){
				trans_bitmap[item] = bitmap_init(bitmap_length);
			}
			//push the tranid into the bitmap
			bitmap_push(trans_bitmap[item], tranid);
		}
		tranid ++;
	}
	//built l1
	list<Itemset>* L1_itemset = new list<Itemset>();
	for(int bcnt=0; bcnt < sz_of_item; bcnt++){
		if( trans_bitmap[bcnt] != NULL){
			Itemset a;
			Itemset_init(a, 1);
			a.items[0] = bcnt;
			a.bitmap = trans_bitmap[bcnt];
			L1_itemset->push_back(a);
		}
	}
	Ln_itemsets.push_back(L1_itemset);
}
//show the content of bitmap
void dumpBitmap(const char* filename){
	uint bmap;
	FILE* pFile = fopen(filename, "wb");
	if (pFile==NULL) perror ("Error opening file");
	else{
		for(int bcnt=0; bcnt < sz_of_item; bcnt++){
			if(trans_bitmap[bcnt] == NULL)
				continue;
			fprintf(pFile, "%d:\n", bcnt);
			for(int dbmap=0; dbmap < bitmap_length; dbmap++){
				bmap = (uint)trans_bitmap[bcnt][dbmap];
				//fprintf(pFile, "0x%x ", trans_bitmap[bcnt][dbmap]);
				int step=0;
				while(bmap != 0){
					if((bmap & 0x00000001)){
						fprintf(pFile, "%d ", int_size * dbmap + step );
					}
					bmap = bmap >> 1;
					step++;
				}
			}
			fprintf(pFile, "\n");
		}
		
		fclose(pFile);
	}

}

//check the n-1 items, used in join two candidate frequent itemset
//level: the current level of the inputted Itemset
bool apriori_chk_ln_1(const Itemset &ck, const Itemset &ck_1, int level){
	//do not check the last one

	int* ck_item = ck.items;
	int* ck_1_item = ck_1.items;
	
	int and_cnt=0;
	int tmp;
	for(int andcnt=0; andcnt < level-1; andcnt++){
		
		if(ck_item[andcnt] != ck_1_item[andcnt]){
				return false;
		}			
	}
	return true;
}
//generate the candidate with gpu
bool apriori_generate_candidate(){

	clock_t gencandidate_str = clock();

	//current level cnt = the length of current candidate items
	int curlevel = Ln_itemsets.size();
	int cur_arr_index = curlevel-1;

	int nextlevel = curlevel+1;
	int next_arr_index = nextlevel -1;

	list<Itemset>*Ln = Ln_itemsets.back();
	if(Ln->size() <2){
		return false;
	}
	list<Itemset>*Lnp1 = new list<Itemset>();
	list<Itemset>::const_iterator c_outter_iter = Ln->begin();
	list<Itemset>::const_iterator c_inner_iter;
	
	for(;c_outter_iter != Ln->end(); ++c_outter_iter){
		c_inner_iter = c_outter_iter;
		c_inner_iter ++;
		for(; c_inner_iter != Ln->end(); ++c_inner_iter){
			//check the Ln_1 is equal or not
			if( apriori_chk_ln_1( (*c_inner_iter), (*c_outter_iter), curlevel )){
				//the previous items are equal, then generate a new lnp1
				Itemset lnp1_item;
				Itemset_init(lnp1_item, nextlevel);
				//sort with insc by check the last item
				if(c_inner_iter->items[cur_arr_index] > c_outter_iter->items[cur_arr_index]){
					//push the outter first and push the inner
					
					
					Itemset_pushArr(lnp1_item, (*c_outter_iter), curlevel);
					
					Itemset_push(lnp1_item, next_arr_index, c_inner_iter->items[cur_arr_index]);
					
					
				}
				else{
					//push the inner first and push the outter

					Itemset_pushArr(lnp1_item, (*c_inner_iter), curlevel);
					
					Itemset_push(lnp1_item, next_arr_index, c_outter_iter->items[cur_arr_index]);

				}
				lnp1_item.bitmap = bitmap_intersect(c_inner_iter->bitmap, c_outter_iter->bitmap, bitmap_length);
				if(lnp1_item.bitmap == NULL){
					//in-frequent
					Itemset_free(lnp1_item);
				}
				else
					Lnp1->push_back(lnp1_item);
			}
		}
	}

	Ln = Ln_itemsets.back();
	//clear Ln's bitmap and items
	list<Itemset>::iterator free_ln_iter = Ln->begin();
	for(; free_ln_iter != Ln->end(); ++free_ln_iter){
		delete[] free_ln_iter->bitmap;
		delete[] free_ln_iter->items;
	}

	//push Lnp1 into Ln_itemsets
	Ln_itemsets.push_back(Lnp1);

	clock_t gencandidate_end = clock();
	fprintf(stderr, "generate candidate cost:%f\n", (double)(gencandidate_end-gencandidate_str)/ (double)CLOCKS_PER_SEC);
	return true;
}
void apriori_prune(){
	clock_t prune_str = clock();
	uint bmap;
	int freq;
	list<Itemset>*Ln = Ln_itemsets.back();
	list<Itemset>::iterator lniter = Ln->begin();
	fprintf(stderr, "prunce no. %d\n",Ln->size());
	while( lniter != Ln->end()){
		freq = 0;
		//count the frequence from the bitmap
		for(int bcnt=0; bcnt < bitmap_length; bcnt++){
			bmap = (uint)lniter->bitmap[bcnt];
			if(bmap == 0)
				continue;
			while(bmap != 0){
				if((bmap & 0x00000001)){
					freq++;
				}
				bmap = bmap >> 1;
			}
			if(freq >= freq_for_minsup)
				break;
		}
		if(freq < freq_for_minsup){
			//delete this
			delete[] lniter->bitmap;
			delete[] lniter->items;
			lniter = Ln->erase(lniter);
		}
		else
			++lniter;
	}
	

	clock_t prune_end = clock();
	fprintf(stderr, "prune cost:%f\n", (double)(prune_end-prune_str)/ (double)CLOCKS_PER_SEC);
}


void apriori_run(){

	cerr << "level:1 "<< Ln_itemsets.back()->size() << endl;
	int level = 2;
	while(apriori_generate_candidate()){
		cerr << "level:" << level << " " << Ln_itemsets.back()->size() << endl;
		level++;
	}
}

int main(int argc, char* argv[]){
	clock_t str = clock();
	//const char* dataset = "AsciiT10.I8.D100K-900.txt";
	const char* dataset = argv[1];
	minsup = atoi(argv[2]);

	loadData(dataset);	
	//dumpTrans("output.tmp");
	convertData2Bitmap();
	//dumpBitmap("output.bitmap.tmp");
	fprintf(stderr, "minsup:%d\n", freq_for_minsup);
	apriori_run();
	freeData();
	
	clock_t end = clock();
	fprintf(stderr, "cost:%f\n", (double)(end-str)/ (double)CLOCKS_PER_SEC);
	
	return 0;
}


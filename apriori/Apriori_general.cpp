/* author: hsinhoyeh yhh92u@gmail.com
 * goal : this cpp implement the basic concept of apriori algorithm
 * */

#include <iostream>
#include <vector>
#include <list>
#include <map>
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
//the max number of items in each frequent itemset
#define sz_of_itemset 32
//the minisup
//#define minsup 0.004
double minsup = 0.01;

typedef unsigned int uint;
typedef struct _itemset{
//change the structure of items into array
	int* items;
	int freq;
} Itemset;

enum prune_type{
	prune_tran,
	prune_itemset,
	insert_itemset,
	insert_nothing
};

list< int* > trans;
list< list<Itemset>* > Ln_itemsets;

//store the L1 frequence
int * L1;
int size_of_trans = -1;
int freq_for_minsup=-1;

void apriori_prune(map<uint, list< Itemset*> > &hashmap);
void apriori_run();
void apriori_subset(int* trans, int* itemset, int level, int nth, int str, map<uint, list< Itemset*> > &hashmap);
bool apriori_generate_candidate();
uint hashvalue(int* itemset, int level);
bool apriori_chk_ln_1(const int* ck, const int* ck_1, int level);

//release the itemset
void Itemset_free(Itemset& a){
	if(a.items != NULL)
		delete[] a.items;
}
//initialize the itemset
void Itemset_init(Itemset& a, int level){
	a.items = new int[level];
	a.freq = 0;
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

		//built l1
		list<Itemset>* L1_itemset = new list<Itemset>();
		for(int bcnt=0; bcnt < sz_of_item; bcnt++){
			if( L1[bcnt] >= freq_for_minsup){
				Itemset a;
				Itemset_init(a, 1);
				a.items[0] = bcnt;
				L1_itemset->push_back(a);
			}
		}
		Ln_itemsets.push_back(L1_itemset);

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

//check the n-1 items, used in join two candidate frequent itemset
//level: the current level of the inputted Itemset
bool apriori_chk_ln_1(const int* ck_item, const int* ck_1_item, int level){
	//do not check the last one
	
	int and_cnt=0;
	for(int andcnt=0; andcnt < level-1; andcnt++){
		if(ck_item[andcnt] != ck_1_item[andcnt]){
			return false;
		}
	}
	return true;
}
//generate the hashvalue with specified itemset
uint hashvalue(int* itemset, int level){
	uint value1=0;
	uint value2=0xffffffff;
	int scale = 1;
	for(int lcnt=0; lcnt< level; lcnt++){
		value1 = value1 + (scale * itemset[lcnt]);
		value2 = value2 - (scale * itemset[lcnt]);
		scale *= 10;
		
	}
	return ( value1 * value2 );
}
//generate candidate
bool apriori_generate_candidate(){

	clock_t gencandidate_str = clock();

	map<uint, list< Itemset* > > hashmap;

	//current level cnt = the length of current candidate items
	int curlevel = Ln_itemsets.size();
	int cur_arr_index = curlevel-1;

	int nextlevel = curlevel+1;
	int next_arr_index = nextlevel -1;

	list<Itemset>*Ln = Ln_itemsets.back();
	//if the number of frequent itemset small than 2, done
	if(Ln->size() <2){		
		return false;
	}
	//else, generate the Cn+1 
	list<Itemset>*Lnp1 = new list<Itemset>();
	list<Itemset>::const_iterator c_outter_iter = Ln->begin();
	list<Itemset>::const_iterator c_inner_iter;
	
	for(;c_outter_iter != Ln->end(); ++c_outter_iter){
		c_inner_iter = c_outter_iter;
		c_inner_iter ++;
		for(; c_inner_iter != Ln->end(); ++c_inner_iter){
			//check the Ln_1 is equal or not
			if( apriori_chk_ln_1( c_inner_iter->items, c_outter_iter->items , curlevel )){
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
				Lnp1->push_back( lnp1_item );
				
				//generate the hash value and push into the hash table
				uint hash_value = hashvalue(lnp1_item.items, nextlevel);
				map<uint, list<Itemset*> >::iterator fiter;
				if( (fiter = hashmap.find(hash_value)) == hashmap.end()){
					list< Itemset* > a;
					a.push_back( &Lnp1->back() );
					hashmap.insert(pair<uint, list< Itemset* >>( hash_value, a ));
				}
				else{
					fiter->second.push_back( &Lnp1->back() );
				}
				
			}
		}
	}

	Ln = Ln_itemsets.back();
	//clear Ln's content
	list<Itemset>::iterator free_ln_iter = Ln->begin();
	for(; free_ln_iter != Ln->end(); ++free_ln_iter){
		delete[] free_ln_iter->items;
	}

	//push Lnp1 into Ln_itemsets
	Ln_itemsets.push_back(Lnp1);

	clock_t gencandidate_end = clock();
	
	fprintf(stderr, "generate candidate cost:%f\n", (double)(gencandidate_end-gencandidate_str)/ (double)CLOCKS_PER_SEC);
	
	//after generate the candidate, then do the prune
	apriori_prune(hashmap);
	
	return true;
}
void apriori_subset(int* trans, int* itemset, int level, int nth, int str, map<uint, list< Itemset*> > &hashmap){
	
	if(str == trans[0] && nth != level){
		//fprintf(stderr, "str == trans[0] && nth != level\n");
		return;
	}
	else if(nth == level){
/*
		fprintf(stderr, "----\n");
		for(int lcnt=0; lcnt< level; lcnt++){
			fprintf(stderr, "%d ", itemset[lcnt]);
		}
		fprintf(stderr, "\n");
*/
		uint hash_value = hashvalue(itemset, level);
		//hashmap
		map<uint, list< Itemset* > >::iterator fiter;
		//increase the frequence with the same hashvalue and same item content
		if( (fiter=hashmap.find(hash_value)) != hashmap.end()){
			list< Itemset* >::iterator siter = fiter->second.begin();
			list< Itemset* >::iterator eiter = fiter->second.end();
			for(; siter != eiter; ++siter){
				if(apriori_chk_ln_1( (*siter)->items, itemset , level+1 )){
					(*siter)->freq++;
				}
			}
		}
		
		return ;
	}
	else{
		for(int tcnt=str; tcnt< trans[0]; tcnt++){
			//pick this one, increase the freq
			
			itemset[nth] = trans[tcnt];
			apriori_subset(trans, itemset, level, nth+1, tcnt+1, hashmap);			
		}
	}
}
//count the itemset's frequence by scaning the whole database
void apriori_prune(map<uint, list< Itemset*> > &hashmap){

	clock_t prune_str = clock();

	list<Itemset>*Ln = Ln_itemsets.back();
	int level = Ln_itemsets.size();

	int* tranItemset = new int[level];
	for(int ini=0; ini<level; ini++){
		tranItemset[ini] = -1;
	}
	prune_type status = insert_nothing;
	//transaction iterator
	list< int* >::iterator titer = trans.begin();
	while(titer != trans.end()){
		status = insert_nothing;

		//scan current transaction
		if( (*titer)[0] <= level ){
			status = prune_tran;
		}
		else{
			apriori_subset((*titer), tranItemset, level, 0, 1, hashmap);
		}


		if(prune_tran == status){
			titer = trans.erase(titer);
		}
		else{
			++titer;
		}

	}
	list<Itemset> prune_donelist;
	//prune
	map<uint, list< Itemset* > >::iterator scan_4_prune_iter = hashmap.begin();
	list< Itemset* >::iterator slistiter;
	list< Itemset* >::iterator elistiter;
	for(; scan_4_prune_iter != hashmap.end(); ++scan_4_prune_iter){
		slistiter = scan_4_prune_iter->second.begin();
		elistiter = scan_4_prune_iter->second.end();
		for(; slistiter != elistiter; ++slistiter){
			if((*slistiter)->freq< freq_for_minsup){
				delete[] (*slistiter)->items;
			}
			else{
				prune_donelist.push_back( (*(*slistiter)) );
			}
		}
	}
	delete[] tranItemset;
	clock_t prune_end = clock();
	
	Ln = Ln_itemsets.back();
	Ln->clear();
	Ln->insert(Ln->end(), prune_donelist.begin(), prune_donelist.end());
	
	fprintf(stderr, "prune cost:%f\n", (double)(prune_end-prune_str)/ (double)CLOCKS_PER_SEC);
}


void apriori_run(){

	cerr << "level:1 "<< Ln_itemsets.back()->size() << endl;
	int level = 2;
	while(true){
		if(apriori_generate_candidate()){
			cerr << "level:" << level << " " << Ln_itemsets.back()->size() << endl;
		}
		else
			break;
		
		level++;
	}
}

int main(int argc, char* argv[]){
	clock_t str = clock();
	//const char* dataset = "AsciiT10.I8.D100K-900.txt";
	const char* dataset = argv[1];
	minsup = atoi(argv[2]);

	loadData(dataset);	
	fprintf(stderr, "minsup:%d\n", freq_for_minsup);
	apriori_run();
	freeData();
	clock_t end = clock();
	fprintf(stderr, "cost:%f\n", (double)(end-str)/ (double)CLOCKS_PER_SEC);
	
	return 0;
}


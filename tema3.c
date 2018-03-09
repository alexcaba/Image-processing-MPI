#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//aplica filtru pentru un pixel din bloc (i,j)
int applyFilterForPixel(int filterTag, int*recv, int height, int width, int i, int j){
	int edgeVec[10] = {1, 1, 0, -1, 2, 0, -2, 1, 0, -1};  //edge detection
	int meanVec[10] = {1, -1, -1, -1, -1, 9, -1, -1, -1, -1};  //mean removal
	int rez;

	if(filterTag == 1){ //filtrul 1
		rez = recv[(i - 1) * width + j - 1] * edgeVec[1] +
			recv[(i - 1) * width + j] * edgeVec[2] +
			recv[(i - 1) * width + j + 1] * edgeVec[3] +
			recv[i * width + j - 1] * edgeVec[4] +
			recv[i * width + j] * edgeVec[5] +
			recv[i * width + j + 1] * edgeVec[6] +
			recv[(i + 1) * width + j - 1] * edgeVec[7] +
			recv[(i + 1) * width + j] * edgeVec[8] +
			recv[(i + 1) * width + j + 1] * edgeVec[9];
		rez += 127;
	}
	else{  //filtrul 2
		rez = recv[(i - 1) * width + j - 1] * meanVec[1] +
			recv[(i - 1) * width + j] * meanVec[2] +
			recv[(i - 1) * width + j + 1] * meanVec[3] +
			recv[i * width + j - 1] * meanVec[4] +
			recv[i * width + j] * meanVec[5] +
			recv[i * width + j + 1] * meanVec[6] +
			recv[(i + 1) * width + j - 1] * meanVec[7] +
			recv[(i + 1) * width + j] * meanVec[8] +
			recv[(i + 1) * width + j + 1] * meanVec[9];
	}

	if (rez > 255){
		return 255;
	}
	if (rez < 0){
		return 0;
	}
	return rez;
}

//aplica filtru pentru un bloc de dimensiuni height x width in functie de filterTag
void applyFilter(int filterTag, int* recv, int height, int width, int* result){
	for (int i = 1; i < height - 1; i++){
		for (int j = 1; j < width - 1; j++){
			result[i * width + j] = applyFilterForPixel(filterTag, recv, height, width, i, j);
		}
	}
}

int main(int argc, char* argv[]){

	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	FILE* top_file = fopen(argv[1], "r");
	int aux, nrVecini = 0;
	int vecini[size];
	
	//citire topologie
	//nodul k sare pana la linia care ii corespunde
	char* line = malloc (4 * size);
	size_t length = 4 * size;

	int lineNo = 0;
	while (lineNo < rank) {
		getline(&line, &length, top_file);
		lineNo++;
	}
	
	fscanf(top_file, "%d:", &aux);
	
	//citire vecini pt nodul curent	
	char temp;
	do{
		fscanf(top_file, "%d%c", &aux, &temp); 
		vecini[nrVecini++] = aux;
	} while(temp!= '\n');

	//determinare topologie
	int i, parent = -1;
	int j, k;
	int number = 10;
	int frunza = 1;

	if(rank == 0){
		//radacina trimite mesaj catre toti vecinii
		for(i = 0; i < nrVecini; i++){
			MPI_Send(&number, 1, MPI_INT, vecini[i], 0, MPI_COMM_WORLD);
			frunza = 0;
		}

		//astept raspuns de la fiecare vecin caruia i-am trimis mesaj
		for(i = 0; i < nrVecini; i++){
			MPI_Recv(&number, 1, MPI_INT, vecini[i], 0, MPI_COMM_WORLD, &status);
		}
	}
	else{
		//astept mesaj de la parinte si il retin
		MPI_Recv(&number, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		if(parent == -1)
			parent = status.MPI_SOURCE;

		//trimit mesajul mai departe
		for(i = 0; i < nrVecini; i++)
			if(vecini[i] != parent){
				MPI_Send(&number, 1, MPI_INT, vecini[i], 0, MPI_COMM_WORLD);
				frunza = 0;
			}

		//astept raspuns de la vecinii carora le-am trimis mesaj
		for(i = 0; i < nrVecini; i++)
			if(vecini[i] != parent){	
				MPI_Recv(&number, 1, MPI_INT, vecini[i], 0, MPI_COMM_WORLD, &status);
			}
		
		//dau raspuns inapoi parintelui
		MPI_Send(&number, 1, MPI_INT, parent, 0, MPI_COMM_WORLD);

		//elimin parintele din vectorul vecini si astfel raman doar cu copiii acestuia
		for(i = 0; i < nrVecini; i++)
			if(vecini[i] == parent) break;
		for(; i < nrVecini; i++)
			vecini[i] = vecini[i + 1];
		nrVecini--;

		//printf("Procesul [%d] are parintele [%d], frunza = %d.\n", rank, parent, frunza);
		//for(i=0;i < nrVecini; i++)
		//	printf("%d ", vecini[i]);
		//printf("\n");
	}

	//procesarea imaginilor
	int nrImg, width, height, maxGV, filterTag = 0;
	int liniiProcesate = 0;
	int* statistica = calloc (size, sizeof(int));
	int* recvStatis = calloc (size, sizeof(int));
	int vect[2];   //folosit pt transmiterea dim blocurilor de procesat
	int vectRec[2];

	if(rank == 0){
		//deschidere fisier imagini.in
		FILE* img_file = fopen(argv[2], "r");
		fscanf(img_file, "%d", &nrImg);
		char filter[256], image1[256], image2[256];

		//retin nume si filtru pt fiecare imagine de procesat
		for(i = 0; i < nrImg; i++){
			fscanf(img_file, "%s", filter);
			fscanf(img_file, "%s", image1);
			fscanf(img_file, "%s", image2);
			
			if(strcmp(filter, "sobel") == 0)
				filterTag = 1;
			if(strcmp(filter, "mean_removal") == 0)
				filterTag = 2;

			//deschidere imagine
			FILE* crtImg = fopen(image1, "r");
			char line1[10];
			char line2[100];
			fgets(line1, 10, crtImg);
			fgets(line2, 100, crtImg);

			fscanf(crtImg, "%d%d%d", &width, &height, &maxGV);

			//citire pixeli
			int* image = calloc((height + 2) * (width + 2), sizeof(int));
			for(j = 1; j <= height; j++)
				for(k = 1; k <= width; k++)
					fscanf(crtImg, "%d", &image[j * (width + 2) + k]);

			//impartire linii de catre radacina
			if(nrVecini > 0){
				int nrVeciniActivi = nrVecini;
				int liniiVecin = height / nrVecini;
				int liniiSend;
				int rest = height % nrVecini;
				int liniaCurenta = 1;
				int primaLinieVecin[size];

				for(j = 0; j < nrVecini; j++){  //trimitere blocuri de prelucrat catre fiecare vecin
					if(rest && (j + 1 == nrVecini))
						liniiSend = liniiVecin + rest;
					else
						liniiSend = liniiVecin;


					vect[0] = liniiSend + 2;
					vect[1] = width + 2;

					if(liniiSend){  //construire bloc de prelucrat
						primaLinieVecin[vecini[j]] = liniaCurenta - 1;

						int* trimite = calloc(vect[0] * vect[1], sizeof(int));
						for(k = (liniaCurenta - 1) * (width + 2); k < (liniaCurenta + liniiSend + 1) * (width + 2); k++)
							trimite[k - (liniaCurenta - 1) * (width + 2)] = image[k];

						MPI_Send(vect, 2, MPI_INT, vecini[j], 0, MPI_COMM_WORLD);
						MPI_Send(trimite, vect[0] * vect[1], MPI_INT, vecini[j], filterTag, MPI_COMM_WORLD);

						liniaCurenta = liniaCurenta + liniiSend;
					}
					else{ //daca nu raman linii de prelucrat anunt fiul
						vect[0] = -1;
						MPI_Send(vect, 2, MPI_INT, vecini[j], 0, MPI_COMM_WORLD);
						nrVeciniActivi--;
					}
				}

				//astept rezultatele de la vecini
				int* rezultat = calloc((height + 2) * (width + 2), sizeof(int));
				for(j = 0; j < nrVeciniActivi; j++){
					MPI_Recv(vect, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					int* primit = calloc (vect[0] * vect[1], sizeof(int));
					MPI_Recv(primit, vect[0] * vect[1], MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				
					//prima linie asignata vecinilor... folosita pt reconstructia imaginii
					liniaCurenta = primaLinieVecin[status.MPI_SOURCE] + 1;
					int l = width + 2;
					for (k = liniaCurenta * (width+2); k < (liniaCurenta+vect[0] - 2) * vect[1]; k++){
						rezultat[k] = primit[l];
						l++;
					}
				}

				//scriere in fisier... poza prelucrata
				FILE* img_file2 = fopen (image2, "w");
				fprintf(img_file2, "%s", line1);
				fprintf(img_file2, "%s", line2);
				fprintf(img_file2, "%d %d\n%d\n", width, height, maxGV);
				for (j = 1; j <= height; j++){
					for (k = 1; k <= width; k++){
						fprintf(img_file2, "%d\n", rezultat[j * (width + 2) + k]);
					}
				}
				fclose(img_file2);
			}
		}

		//trimitere mesaj cu tag final dupa prelucrarea pozelor
		for (i = 0; i < nrVecini; i++)
			MPI_Send(vect, 2, MPI_INT, vecini[i], 5, MPI_COMM_WORLD);

		//astept statistica de la fiecare vecin si asamblez rezultatele
		for(i = 0; i < nrVecini; i++){
			MPI_Recv(recvStatis, size, MPI_INT, vecini[i], MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			int m;
			for(m = 0; m < size; m++)
				if(recvStatis[m] > statistica[m])
					statistica[m] = recvStatis[m];
		}

		//scriere statistica in fisier
		FILE* statFile = fopen(argv[3], "w");
		for(int x = 0; x < size; x++)
			fprintf(statFile, "%d: %d\n", x, statistica[x]);
		fclose(statFile);
	}

	if(rank != 0){

		while(1){ //cat timp nu am primit tag-ul de incheiere
			int* vectRec = calloc (2, sizeof(int)); 
			MPI_Recv(vectRec, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			//tag final
			if(status.MPI_TAG == 5){
				//trimit mesajul de final mai departe
				for(i = 0; i < nrVecini; i++)
					MPI_Send(vectRec, 2, MPI_INT, vecini[i], 5, MPI_COMM_WORLD);

				//daca sunt frunza trimit statisica
				if(frunza == 1){
					statistica[rank] = liniiProcesate;
					MPI_Send(statistica, size, MPI_INT, parent, 0, MPI_COMM_WORLD);
				}
				else{ //altfel astept statistica de la fiecare fiu, asamblez rezultatele si trimit mai departe
					for(i = 0; i < nrVecini; i++){
						MPI_Recv(recvStatis, size, MPI_INT, vecini[i], MPI_ANY_TAG, MPI_COMM_WORLD, &status);
						int m;
						for(m = 0; m < size; m++)
							if(recvStatis[m] > statistica[m])
								statistica[m] = recvStatis[m];
					}
					MPI_Send(statistica, size, MPI_INT, parent, 0, MPI_COMM_WORLD);
				}
				break;
			}

			height = vectRec[0];
			width = vectRec[1];
			if(height > 0){ //daca am primit linii de procesat

				int* primit = calloc(height * width, sizeof(int));

				MPI_Recv(primit, height * width, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
				filterTag = status.MPI_TAG;

				if(nrVecini > 0){ //nodul nu este frunza, impart mai departe
					int nrVeciniActivi = nrVecini;
					int liniiVecin = (height - 2) / nrVecini;
					int liniiSend;
					int rest = (height - 2) % nrVecini;
					int liniaCurenta = 1;
					int primaLinieVecin[size];

					for (j = 0; j < nrVecini; j++) {  //trimit la fiecare vecin
							if (rest && (j + 1 == nrVecini)) {
								liniiSend = liniiVecin + rest;
							}
							else {
								liniiSend = liniiVecin;
							}
							
							//construiesc bloc de trimis		
							if (liniiSend > 0) {

								vect[0] = liniiSend + 2;
								vect[1] = width;
								
								int* trimite = calloc(vect[0] * vect[1], sizeof(int));
								int l = 0;
								primaLinieVecin[vecini[j]] = liniaCurenta - 1;
								for (k = (liniaCurenta - 1) * width; k < (liniaCurenta + liniiSend + 1) * width; k++){
									trimite[l] = primit[k];
									l++;
								}

								MPI_Send(vect, 2, MPI_INT, vecini[j], 0, MPI_COMM_WORLD);
								MPI_Send(trimite, vect[0] * vect[1], MPI_INT, vecini[j], filterTag, MPI_COMM_WORLD);
								liniaCurenta += liniiSend;
							}
							else {  //daca nu exista nimic de prelucrat anunt vecinii
								vect[0] = -1;
								MPI_Send(vect, 2, MPI_INT, vecini[j], 1, MPI_COMM_WORLD);
								nrVeciniActivi--;
							}
					}

					int* rezultat = calloc (height * width, sizeof(int));

					//astept rezultatele de la vecinii care au primit linii de procesat si le combin
					for (j = 0; j<nrVeciniActivi; j++){
						MPI_Recv(vectRec, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
						primit = calloc(vectRec[0] * vectRec[1], sizeof(int));
						MPI_Recv(primit, vectRec[0] * vectRec[1], MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

						liniaCurenta = primaLinieVecin[status.MPI_SOURCE] + 1;
						int l = width;
						
						for (k = liniaCurenta * width; k<(liniaCurenta+vectRec[0] - 2)*vectRec[1]; k++){
							rezultat[k] = primit[l];
							l++;	
						}
					}

					//trimit inapoi ce au procesat frunzele
					vectRec[0] = height;
					vectRec[1] = width;
					MPI_Send(vectRec, 2, MPI_INT, parent, 0, MPI_COMM_WORLD);
					MPI_Send(rezultat, height * width, MPI_INT, parent, 0, MPI_COMM_WORLD);
				}
				else{ //nodul este frunza si aplica filtrul pe blocul primit

					int* prelucrare = calloc(height * width, sizeof(int));
					memcpy(prelucrare, primit, height * width * sizeof(int));

					applyFilter(filterTag, primit, height, width, prelucrare);
					liniiProcesate += height - 2;  //numara liniile procesate

					MPI_Send(vectRec, 2, MPI_INT, parent, 0, MPI_COMM_WORLD);
					MPI_Send(prelucrare, height * width, MPI_INT, parent, 0, MPI_COMM_WORLD);
				}
			}
			else{ //nu are linii de prelucrat
				vectRec[0] = -1;
				for (i = 0; i < nrVecini; i++){
					if(vecini[i] != parent)
						MPI_Send(vectRec, 2, MPI_INT, vecini[i], 0, MPI_COMM_WORLD);
				}
			}
		}
	}

	MPI_Finalize();
	return 0;
}
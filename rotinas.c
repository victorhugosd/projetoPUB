#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int levels_in_recip[10001];

typedef struct {
    double x, y;
} point;

struct Point_node {
    point p;
    struct Point_node *next;
};

typedef struct Point_node point_node;

typedef struct Index_node {
    int i;
    struct Index_node *next;
} index_node;

typedef index_node* level_list;

typedef struct {
    double r, h;
    point center; //ponto no plano da base do recipiente em que o cilindro será empacotado
} cylinder;

point_node*newPoint_list () {
    return NULL;
}

//insere ponto no início de uma lista
point_node *insert_point_in_list (point_node *list_points, point P) {
    point_node *pointer = (point_node*)malloc(sizeof(point_node));
    if(pointer!=NULL) {
        pointer->p = P;
        pointer -> next = list_points;
        return pointer;
    }
    return list_points;
}

void print_points_in_list (point_node *list_points) {
    point_node *pointer;
    for (pointer = list_points; pointer!=NULL; pointer = pointer->next) {
        printf("(%lf , %lf)\n", pointer->p.x, pointer->p.y);
    }
    return;
}

index_node *newIndex_list () {
    return NULL;
}

index_node *insert_index_in_list (index_node *list_indices, int i) {
    index_node *pointer = (index_node*)malloc(sizeof(index_node));
    pointer->i = i;
    if(pointer!=NULL) {
        pointer->i = i;
        pointer -> next = list_indices;
        return pointer;
    }
    return list_indices;
}

void print_indices_in_list (index_node *list_indices) {
    index_node *pointer;
    for (pointer = list_indices; pointer!=NULL; pointer = pointer->next) {
        printf("%d ", pointer->i);
    }
    printf("\n");
    return;
}

double max (double a, double b) {
    return (a > b ? a : b);
}

double min (double a, double b) {
    return (a < b ? a : b);
}

double abs_value (double a) {
    return (a<0 ? -a : a);
}

//verifica se dois valores são iguais dentro de uma faixa de tolerância
int equals (double a, double b, double L, double W) {
    double d = max(L, W);
    if (abs_value(a-b)<max(1e-5, d*(1e-8))) return 1;
    return 0;
}

//verifica se um valor (primeiro parrÂmetro) é menor que outro além de uma faixa tolerável de acordo com as dimensões do recipiente
int less_than (double a, double b, double L, double W) {
    double d = max(L, W);
    if (a < b - max(1e-5, d*(1e-8))) return 1;
    return 0;
}

//verifica se um valor (primeiro parrÂmetro) é maior que outro além de uma faixa tolerável
int greater_than (double a, double b, double L, double W) {
    double d = max(L, W);
    if (a > b + max(1e-5, d*(1e-8))) return 1;
    return 0;
}

//retorna a distância entre dois pontos
double distance (point p1, point p2) {
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

//verifica se o cilindro está dentro da região de contenção do recipiente
int in_container (cylinder c, point p, double L, double W) {
	if (less_than(p.x, c.r, L, W) || greater_than(p.x, L-c.r, L, W) || less_than(p.y, c.r, L, W) || greater_than(p.y, W-c.r, L, W)) {
		return 0;
	}
	return 1;
}

//verifica se há sobreposição do cilindro c com os cilindros empacotados anteriormente
int overlap (point p, cylinder *list_itens, int index, index_node *level, double L, double W) {
	int i;
	index_node *pointer = level;
	while (pointer) {
        i = pointer->i;
        if (less_than(distance(p, list_itens[i].center), list_itens[index].r + list_itens[i].r, L, W)) {
            //printf("Nao, pois se empacotado nessa posicao, %d se sobrepoe ao item %d\n", index, i);
            return 1;
        }
        pointer = pointer->next;
	}
	return 0;
}

//verifica se um ponto p é fáctível para empacotar a base do cilindro de índice index no nível level
int is_feasible_point (point p, cylinder* list_itens, int index, index_node *level, double L, double W) {
    if (in_container (list_itens[index], p, L, W)) {
        if (!overlap(p, list_itens, index, level, L, W)) return 1;
        else return 0;
    }
    return 0;
}

//verifica se o ponto p2 está à esquerda, além de uma tolerância, do ponto p1
//se ambos tiverem mesma coordenada x,  utiliza-se a menor coordenada y como critério de desempate
int more_to_left (point p1, point p2, double L, double W) {
        if (less_than(p2.x, p1.x, L, W)) {
            return 1;
        }
        else if (equals(p2.x, p1.x, L, W)) {
            if (less_than(p2.y, p1.y, L, W)) {
                return 1;
            }
        }
        return 0;
}

//encontra a interseção das linhas verticais da fronteira da região de contenção com os a fronteira da região de sobreposição com itens já alocados no nível
int find_inter_vertical_line_circle (point *p1, point *p2, point c, double x0, double r1, double r, double L, double W) {
    double D = (r1 + r)*(r1 + r) - (x0 - c.x)*(x0 - c.x);
    if (equals(D, 0.0, L, W)) {
        p1->x = x0;
        p1->y = c.y;
        return 1;
    }
    else if (greater_than(D, 0.0, L, W)) {
        p1->x = p2->x = x0;
        p1->y = c.y + sqrt(D);
        p2->y = c.y - sqrt(D);
        return 2;
    }
    else return 0;
}

//encontra a interseção das linhas horizontais da fronteira da região de contenção com os a fronteira da região de sobreposição com os itens já alocados no nível
int find_inter_horizontal_line_circle(point *p1, point *p2, point c, double y0, double r1, double r, double L, double W) {
    double D = (r1 + r)*(r1 + r) - (y0 - c.y)*(y0 - c.y);
    if (equals(D, 0.0, L, W)) {
        p1->y = y0;
        p1->x = c.x;
        return 1;
    }
    else if (greater_than(D, 0.0, L, W)) {
        p1->y = p2->y = y0;
        p1->x = c.x + sqrt(D);
        p2->x = c.x - sqrt(D);
        return 2;
    }
    else return 0;
}

//encontra os pontos de interseção das fronteiras das regiões de sobreposição de dois dos itens já alocados no nível
int find_inter_two_circles (point *p1, point *p2, point c1, point c2, double r1, double r2, double r, double L, double W) {
    if (greater_than(distance(c1, c2), r1 + r2 + 2*r, L, W)) return 0;
    else if (equals(distance(c1, c2), r1 + r2 + 2*r, L, W)) {
        p1->x = c1.x + (c2.x - c1.x)*(r1 + r)/distance(c1, c2);
        p1->y = c1.y + (c2.y - c1.y)*(r1 + r)/distance(c1, c2);
        return 1;
    }
    else if (less_than(distance(c1, c2), r1 + r2 + 2*r, L, W)) {
        if (equals(c1.y, c2.y, L, W)) {
            if (equals(c1.x, c2.x, L, W)) {
                //printf("Circulos empacotados no mesmo ponto\n");
                return -1;
            }
            else {
                p1->x = p2->x = (c1.x + c2.x + ((r1+r)*(r1+r) - (r2+r)*(r2+r))/(c2.x - c1.x))/2;
                double D = (r1+r)*(r1+r) - (p1->x - c1.x)*(p1->x - c1.x);
                if (D >= 0) {
                    p1->y = c1.y + sqrt(D);
                    p2->y = c2.y + sqrt(D);
                    return 2;
                }
                else {
                    //printf("Delta negativo, algo deu errado\n");
                    return -1;
                }
            }
        }
        else {
            double A = (c1.x - c2.x)/(c2.y - c1.y);
            double B = (c1.y + c2.y + ((r1 + r)*(r1 + r) - (r2 + r)*(r2 + r) + (c2.x)*(c2.x) - (c1.x)*(c1.x))/(c2.y - c1.y))/2;
            //printf("A = %lf B = %lf\n", A, B);

            double a = 1 + A*A;
            //printf("a = %lf\n", a);

            double b = 2*(A*B - c1.x - A*(c1.y));
            //printf("b = %lf", b);

            double c = (c1.x)*(c1.x) + (B - c1.y)*(B - c1.y) - (r1 + r)*(r1 + r);
            //printf("c = %lf", c);

            double D = b*b - 4*a*c;
            //printf("Delta: %lf\n", D);
            if (D >= 0) {
                p1->x = (-b + sqrt(D))/(2*a);
                p1->y = A*(p1->x) + B;

                p2->x = (-b - sqrt(D))/(2*a);
                p2->y = A*(p2->x) + B;
                return 2;
            }
            else {
                //printf("Delta negativo, algum erro\n");
                return -1;
            }
        }
    }
}

//encontra os pontos factíveis de empacotamento segundo a heurística
point *find_most_left (cylinder *list_itens, int index, index_node *levels[], int *num_levels, int *K, double L, double W) {
    int i, j;
    point p, p2, *best_point = (point*)malloc(sizeof(point));
    if (best_point!=NULL) {
        for(*K = 0; *K<*num_levels; (*K)++) {
            int k = *K;
            //printf("*****************************************************************************************8\n");
            //printf("Tentando empacotar o item %d, de raio %lf, no nivel %d\n", index, list_itens[index].r, k);
            //printf("Linha 157:\n");
            //printf("Itens empacotados nesse nivel: ");
            ////print_indices_in_list(levels[k]);

            best_point->x = L;
            best_point->y = W;

            //primeiro, verifica os vértices do polígono de contenção OBS: TESTADO E FUNCIONANDO
            //printf("Veriicando pontos de intersecao das arestas da fronteira da regiao de contencao:\n");
            //printf("-----------------------------------------------------------------\n");

            p.x = list_itens[index].r;
            p.y = W - list_itens[index].r;
            //printf("Ponto calculado: (%lf, %lf) (melhor ate agora: (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
            if (is_feasible_point(p, list_itens, index, levels[k], L, W)) {
                if (more_to_left(*best_point, p, L, W)) {
                    best_point->x = p.x;
                    best_point->y = p.y;
                    //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                }
                //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            }
            p.x = L - list_itens[index].r;
            p.y = W - list_itens[index].r;
            //printf("Ponto calculado: (%lf, %lf) (melhor ate agora: (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
            if (is_feasible_point(p, list_itens, index, levels[k], L, W)) {
                if (more_to_left(*best_point, p, L, W)) {
                    best_point->x = p.x;
                    best_point->y = p.y;
                    //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                }
                //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            }
            p.x = p.y = list_itens[index].r;
            //printf("Ponto calculado: (%lf, %lf) (melhor ate agora: (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
            if (is_feasible_point(p, list_itens, index, levels[k], L, W)) {
                if (more_to_left(*best_point, p, L, W)) {
                    best_point->x = p.x;
                    best_point->y = p.y;
                    //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                }
                //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            }
            p.x = L - list_itens[index].r;
            p.y = list_itens[index].r;
            //printf("Ponto calculado: (%lf, %lf) (melhor ate agora: (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
            if (is_feasible_point(p, list_itens, index, levels[k], L, W)) {
                if (more_to_left(*best_point, p, L, W)) {
                    best_point->x = p.x;
                    best_point->y = p.y;
                    //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                }
                //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
            }

            //depois, a interseção das arestas do polígono de contenção do recipiente com a região de sobreposição com os circulos já empacotados
            //OBS: TAMBÉM FUNCIONADNDO!!!
            //printf("Intersecao das arestas da regiao de contencao com os cilindros ja empacotados:\n");
            //printf("Lembrando que estamos tentando empacotar o item %d no nivel %d, e os itens ja empacotaos nesse nivel sao: ", index, k);
            //print_indices_in_list(levels[k]);
            //printf("-----------------------------------------------\n");

            index_node *pointer1 = NULL;
            index_node *pointer2 = NULL;
            //int cont = 0;
            for (pointer1 = levels[k]; pointer1!=NULL; pointer1 = pointer1->next) {
                i = pointer1->i;
                double R = list_itens[index].r+list_itens[i].r;
                //sobre a reta vertical à esquerda
                //printf(" Checando a reta vertical a esquerda com item %d:\n\n", i);
                int d = find_inter_vertical_line_circle(&p, &p2, list_itens[i].center, list_itens[index].r, list_itens[i].r, list_itens[index].r, L, W);
                if(d == 1) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
                else if (d == 2) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p2, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p2, L, W)) {
                            best_point->x = p2.x;
                            best_point->y = p2.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }

                //vertical à direita
                //printf("Agora, aresta vertical a direita com o item %d:\n\n", i);

                d = find_inter_vertical_line_circle(&p, &p2, list_itens[i].center, L - list_itens[index].r, list_itens[i].r, list_itens[index].r, L, W);;
                if (d == 1) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
                else if (d == 2) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p2, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p2, L, W)) {
                            best_point->x = p2.x;
                            best_point->y = p2.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }

                //horizontal inferior
                //printf("Agora, aresta horizontal inferior com o item %d:\n\n", i);

                d = find_inter_horizontal_line_circle(&p, &p2, list_itens[i].center, list_itens[index].r, list_itens[i].r, list_itens[index].r, L, W);
                if (d == 1) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
                else if (d == 2) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p2, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p2, L, W)) {
                            best_point->x = p2.x;
                            best_point->y = p2.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }

                //horizontal superior
                //printf("Agora, aresta horizontal superior com o item %d:\n\n", i);

                d = find_inter_horizontal_line_circle(&p, &p2, list_itens[i].center, W - list_itens[index].r, list_itens[i].r, list_itens[index].r, L, W);
                if (d == 1) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
                else if (d == 2) {
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p, L, W)) {
                            best_point->x = p.x;
                            best_point->y = p.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                    //printf("Eh possivel empacotar o item %d nesse ponto? ", index);
                    if (is_feasible_point(p2, list_itens, index, levels[k], L, W) && more_to_left(*best_point, p2, L, W)) {
                            best_point->x = p2.x;
                            best_point->y = p2.y;
                            //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);

                        //else printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
            }

            //agora, as interseções dois a dois das regiões de sobreposição dos cilindros já empacotados
            //printf("Verificando as intersecoes dois a dois das regioes de contencao com os cilindros ja empacotados no nivel %d:\n", k);
            //printf("----------------------------------------------------------------\n");

            for (pointer1 = levels[k]; pointer1!=NULL; pointer1 = pointer1->next)
            for (pointer2 = pointer1->next; pointer2!=NULL; pointer2 = pointer2->next) {
                i = pointer1->i;
                j = pointer2->i;
                double x1 = list_itens[i].center.x, x2 = list_itens[j].center.x;
                double y1 = list_itens[i].center.y, y2 = list_itens[j].center.y;
                double r1 = list_itens[i].r+list_itens[index].r, r2 = list_itens[j].r+list_itens[index].r;
                //printf("OS ITENS %d E %d ESTAO EMPACOTADOS NOS PONTOS (%lf, %lf) E (%lf, %lf), RESPECTIVAMENTE\n", i, j, x1, y1, x2, y2);
                int d = find_inter_two_circles(&p, &p2, list_itens[i].center, list_itens[j].center, list_itens[i].r, list_itens[j].r, list_itens[index].r, L, W);
                if (d == 1) {
                    if (more_to_left(*best_point, p, L, W)) {
                        best_point->x = p.x;
                        best_point->y = p.y;
                        //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                    }
                    else {
                        //printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
                else if (d == 2) {
                    if (more_to_left(*best_point, p, L, W) && is_feasible_point(p, list_itens, index, levels[k], L, W)) {
                        best_point->x = p.x;
                        best_point->y = p.y;
                        //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                    }
                    else {
                        //printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }

                    if (more_to_left(*best_point, p2, L, W) && is_feasible_point(p2, list_itens, index, levels[k], L, W)) {
                        best_point->x = p2.x;
                        best_point->y = p2.y;
                        //printf("Este ponto eh melhor que o atual, entao o melhor potno ate agora eh: (%lf, %lf)\n", best_point->x, best_point->y);
                    }
                    else {
                        //printf("Embora seja um ponto factivel, (%lf, %lf) nao eh melhor que o anterior (%lf, %lf)\n", p.x, p.y, best_point->x, best_point->y);
                    }
                }
            }

            if (!equals(best_point->x, L, L, W)) {
                //printf("O item %d foi empacotado no nivel %d no ponto (%lf, %lf)\n\n", index, k, best_point->x, best_point->y);
                return best_point;
            }
            //else printf("O item %d nao pode ser empacotado no nivel %d\n", index, k);
        }
        //printf("Foi criado um novo nivel (nivel %d)\n", *num_levels);
        (*num_levels)++;
        best_point->x = best_point->y = list_itens[index].r;
        //printf("O item %d foi empacotado nesse nivel no ponto (%lf, %lf)\n\n", index, best_point->x, best_point->y);
        return best_point;
    }
    return NULL;
}

//ordena o vetor v de indices sefundo os valores do vetor space (ordem crescente)
void sort_increasing(int v[], double space[], int m, double H) {
    int i, j;
    for (i = 0; i < m - 1; i++)
    for (j = i + 1; j < m; j++) {
        if (less_than(space[v[j]], space[v[i]], H, H)) {
            int temp = v[j];
            v[j] = v[i];
            v[i] = temp;
        }
    }
    return;
}

//ordena o vetor v de indices sefundo os valores do vetor space (ordem decrescente)
void sort_decreasing(int v[], double space[], int m, double H) {
    int i, j;
    for (i = 0; i < m - 1; i++)
    for (j = i + 1; j < m; j++) {
        if (greater_than(space[v[j]], space[v[i]], H, H)) {
            int temp = v[j];
            v[j] = v[i];
            v[i] = temp;
        }
    }
    return;
}

//heurística Best-Fit para alocação dos níveis nos recipientes
level_list *best_fit (double position_level[], double level_height[], int *num_recipients, int num_levels, double H) {
    level_list *ans = (level_list*)malloc(100000*sizeof(level_list));
    int j;
    for (j = 0; j < 100000; j++) {
        ans[j] = newIndex_list();
    }
    *num_recipients = 1;
    double *remaining_space = (double*)malloc(100000*sizeof(double));
    int *less_empty = (int*)malloc(100000*sizeof(int));
    for (j = 0; j < 100000; j++) {
        less_empty[j] = j;
    }
    remaining_space[*num_recipients-1] = H;
    for (j = 0; j < num_levels; j++) {
        //printf("----------------------------------------------\n");
        //printf("Vamos alocar o nivel %d\n", j);
        int i = 0, aloc = 0;
        while (i < *num_recipients && aloc == 0) {
            int k = less_empty[i];
            //printf("Espaco vertical restante no recipiente %d: %lf\n", k, remaining_space[k]);
            //printf("Altura do nivel %d: %lf\n", j, level_height[j]);
            if (!less_than(remaining_space[k], level_height[j], H, H)) {
                aloc = 1;
                //printf("Nivel %d alocado no reipiente %d\n\n", j, k);
                levels_in_recip[k]++;
            }
            else i++;
        }
        if (!aloc) {
            (*num_recipients)++;
            remaining_space[(*num_recipients)-1] = H;
            //printf("O nivel %d não coube em nenhum recipiente existente, foi criado um novo: %d\n", j, *num_recipients-1);
            //printf("Nivel %d alocado no recipiente %d\n\n", j, *num_recipients-1);
            levels_in_recip[*num_recipients-1]++;
        }
        position_level[j] = H - remaining_space[less_empty[i]];
        remaining_space[less_empty[i]] = remaining_space[less_empty[i]] - level_height[j];
        insert_index_in_list(ans[less_empty[i]], j);
        sort_increasing(less_empty, remaining_space, *num_recipients, H);
        //printf("Recipientes em ordem crescente de espaco restante:\n");
        for (i = 0; i < *num_recipients; i++) {
            //printf("Recipiente %d Espaco restante: %lf\n", less_empty[i], remaining_space[less_empty[i]]);
        }
    }
    return ans;
}


//heurística Worst-Fit para alocação dos níveis nos recipientes
level_list *worst_fit (double position_level[], double level_height[], int *num_recipients, int num_levels, double H) {
    level_list *ans = (level_list*)malloc(100000*sizeof(level_list));
    int j;
    for (j = 0; j < 100000; j++) {
        ans[j] = newIndex_list();
    }
    *num_recipients = 1;
    double *remaining_space = (double*)malloc(100000*sizeof(double));
    int *more_empty = (int*)malloc(100000*sizeof(int));
    for (j = 0; j < 100000; j++) {
        more_empty[j] = j;
    }
    remaining_space[*num_recipients-1] = H;
    for (j = 0; j < num_levels; j++) {
        //printf("----------------------------------------------\n");
        //printf("Vamos alocar o nivel %d\n", j);
        int i = 0, aloc = 0;
        while (i < *num_recipients && aloc == 0) {
            int k = more_empty[i];
            //printf("Espaco vertical restante no recipiente %d: %lf\n", k, remaining_space[k]);
            //printf("Altura do nivel %d: %lf\n", j, level_height[j]);
            if (!less_than(remaining_space[k], level_height[j], H, H)) {
                aloc = 1;
                //printf("Nivel %d alocado no reipiente %d\n\n", j, k);
                levels_in_recip[k]++;
            }
            else i++;
        }
        if (!aloc) {
            (*num_recipients)++;
            remaining_space[(*num_recipients)-1] = H;
            //printf("O nivel %d não coube em nenhum recipiente existente, foi criado um novo: %d\n", j, *num_recipients-1);
            //printf("Nivel %d alocado no recipiente %d\n\n", j, *num_recipients-1);
            levels_in_recip[*num_recipients-1]++;
        }
        position_level[j] = H - remaining_space[more_empty[i]];
        remaining_space[more_empty[i]] = remaining_space[more_empty[i]] - level_height[j];
        insert_index_in_list(ans[more_empty[i]], j);
        sort_decreasing(more_empty, remaining_space, *num_recipients, H);
        //printf("Recipientes em ordem decrescente de espaco restante:\n");
        for (i = 0; i < *num_recipients; i++) {
            //printf("Recipiente %d Espaco restante: %lf\n", more_empty[i], remaining_space[more_empty[i]]);
        }
    }
    return ans;
}


//heurística First-Fit para alocação dos níveis nos recipientes
level_list *first_fit (double position_level[], double level_height[], int *num_recipients, int num_levels, double H) {
    level_list *ans = (level_list*)malloc(100000*sizeof(level_list));
    int j;
    for (j = 0; j < 100000; j++) {
        ans[j] = newIndex_list();
    }
    *num_recipients = 1;
    double *remaining_space = (double*)malloc(100000*sizeof(double));
    remaining_space[*num_recipients-1] = H;
    for (j = 0; j < num_levels; j++) {
        //printf("----------------------------------------------\n");
        //printf("Vamos alocar o nivel %d\n", j);
        int i = 0, aloc = 0;
        while (i < *num_recipients && aloc == 0) {
            //printf("Espaco vertical restante no recipiente %d: %lf\n", i, remaining_space[i]);
            //printf("Altura do nivel %d: %lf\n", j, level_height[j]);
            if (!less_than(remaining_space[i], level_height[j], H, H)) {
                aloc = 1;
                //printf("Nivel %d alocado no recipiente %d\n\n", j, i);
                levels_in_recip[i]++;
            }
            else i++;
        }
        if (!aloc) {
            (*num_recipients)++;
            remaining_space[(*num_recipients)-1] = H;
            //printf("O nivel %d não coube em nenhum recipiente existente, foi criado um novo: %d\n", j, *num_recipients-1);
            //printf("Nivel %d alocado no recipiente %d\n\n", j, *num_recipients-1);
            levels_in_recip[*num_recipients-1]++;
        }
        position_level[j] = H - remaining_space[i];
        remaining_space[i] = remaining_space[i] - level_height[j];
        insert_index_in_list(ans[i], j);
    }
    return ans;
}

void print_Tex_levels (cylinder list_itens[], level_list levels[], int num_levels, double L, double W) {
    int j;
    printf("Codigo LaTex:\n\n");
	printf("\\documentclass{article}\n");
    printf("\\usepackage{tikz}\n");
    printf("\\begin{document}\n");
	for (j = 0; j<num_levels; j++) {
        printf("\\begin{tikzpicture}[scale=3]\n");
        printf("\\draw (0,0) -- (4,0) -- (4,%lf) -- (0,%lf) -- cycle;\n", 4*W/L, 4*W/L);
        index_node *pointer = levels[j];
        while (pointer) {
            int i;
            i = pointer->i;
            printf("\\draw (%lf,%lf) circle (%lf) node {$%d (r=%.2lf,x=%.2lf,y=%.2lf)$};\n", 4*list_itens[i].center.x/L, 4*list_itens[i].center.y/L, 4*list_itens[i].r/L, i, list_itens[i].r, list_itens[i].center.x, list_itens[i].center.y);
            pointer = pointer->next;
        }
        printf("\\end{tikzpicture}\n");
        printf("\\newpage\n");
	}
	printf("\\end{document}\n\n");
    return;
}

void print_tab (int a, int b) {
    printf("\nSolucao BL");
    if (a==1) printf("h");
    else if (a==2) printf("r");
    if (b==1) printf("BF:\n\n");
    else if (b==2) printf("FF:\n\n");
    else if (b==3) printf("WF:\n\n");
    return;
}

int main () {

    int i, j, k, n;
    double L, W, H;
    level_list *levels = (level_list*)malloc(1000000*sizeof(level_list));
    for (i = 0; i < 1000000; i++) {
        levels[i] = newIndex_list();
    }
    int num_levels = 1;
	printf("Insira a quantidade de cilindros a serem empacotados:\n");
	scanf("%d\n", &n);
	int *indices = (int*)malloc(n*sizeof(int));
    for (i = 0; i < n; i++) {
        indices[i] = i;
    }
	printf("Insira as dimensões do recipiente:\n");
	printf("Largura: ");
	scanf("%lf\n", &W);
	printf("Comprimento: ");
	scanf("%lf\n", &L);
	printf("Altura: ");
	scanf("%lf\n", &H);
	printf("Insira o tamanho do raio da base e a altura dps itens:\n");
	cylinder *list_itens = (cylinder*)malloc(n*sizeof(cylinder));
	getchar();
	for (i = 0; i < n; i++) {
		(list_itens+i)->center.x = (list_itens+i)->center.y = 0.0;
		scanf("%lf", &list_itens[i].r);
		getchar();
		if (greater_than(2*list_itens[i].r, L, L, W) || greater_than(2*list_itens[i].r, W, L, W) || greater_than(list_itens[i].h, H, H, H)) {
            printf("Um dos itens nao cabe no recipiente vazio\n");
            return 0;
		}
	}
	getchar();

	getchar();
	for (i = 0; i < n; i++) {
		scanf("%lf", &list_itens[i].h);
		getchar();
		if (greater_than(2*list_itens[i].r, L, L, W) || greater_than(2*list_itens[i].r, W, L, W) || greater_than(list_itens[i].h, H, H, H)) {
            printf("Um dos itens nao cabe no recipiente vazio\n");
            return 0;
		}
	}
	getchar();

    unsigned code1, code2;
    line718:
    printf("Como quer ordenar os itens?\n");
    printf("1 - Ordem decrescente de raio\n");
    printf("2 - Ordem decrescente de altura\n");
    printf("3 - Nao ordenar\n");
    scanf("%u", &code);

    for (code1 = 1; code1 < 3; code1++) {
        for (code2 = 1; code2 < 4; code2++) {
            clock_t t = clock();
            if (code1 == 2) {
                printf("Opcao escolhida: Ordem decrescente de raio\n");
                double *radius = (double*)malloc(n*sizeof(double));
                for (i = 0; i < n; i++) {
                    radius[i] = list_itens[i].r;
                }
                sort_decreasing(indices, radius, n, min(L, W));
            }
            else if (code1 == 1) {
                printf("Opcao escolhida: Ordem decrescente de altura\n");
                double *height = (double*)malloc(n*sizeof(double));
                for (i = 0; i < n; i++) {
                    height[i] = list_itens[i].h;
                }
                sort_decreasing(indices, height, n, min(L, W));
            }
            else if (code1 == 3) {
                printf("Opcao escolhida: nao ordenar");
            }
            else {
                printf("Opcao invalida, digite de novo\n");
                goto line718;
            }
        	//Alocando as bases dos itens nos níveis
        	for (i = 0; i < n; i++) {
        	    //printf("Item %d  Raio: %lf\n", i, list_itens[i].r);
        	    int n = num_levels;
        	    point *a = find_most_left(list_itens, indices[i], levels, &num_levels, &k, L, W);
        	    if (a != NULL) {
                    list_itens[indices[i]].center = *a;
                    levels[k] = insert_index_in_list(levels[k], indices[i]);
                }
                else return -1;
        	}
            /*
        	for (j = 0; j<num_levels; j++) {
                printf("Itens empacotados no nivel %d: ", j);
                index_node *pointer = levels[j];
                while (pointer) {
                    printf("%d, ", pointer->i);
                    pointer = pointer->next;
                }
                printf("\n\n");
        	}
            */
        	
            printf("Numero de niveis: %d\n", num_levels);

            print_Tex_levels(list_itens, levels, num_levels, L, W);
            

        	//Agora, alocando os níveis nos recipientes
        	//printf("Quantidade de niveis: %d\n", num_levels);
        	
            double *level_height = (double*)malloc(num_levels*sizeof(double));
        	
            for(i = 0; i < num_levels; i++) {
                //printf("--------\n");
                //printf("Nivel %d\n", i);
                index_node *pointer = levels[i];
                double largest = list_itens[pointer->i].h;
                while (pointer) {
                    int j = pointer->i;
                    //printf("Altura: %lf\n", list_itens[j].h);
                    largest = max(largest, list_itens[j].h);
                    pointer = pointer->next;
                }
                level_height[i] = largest;
                //printf("Altura do nivel %d: %lf\n\n", i, level_height[i]);
        	}

            int num_recipients = 1; //quantidade de recipientes;
            double *position_level = (double*)malloc(num_levels*sizeof(double)); //posicao vertical de cada nível

            printf("************************************************************************************\n");
            printf("************************************************************************************\n");
            printf("Parte 2 da heuristica\n");

            line804:

            printf("Escolha o metodo para alocar os niveis:\n");
            printf("1 - First Fit\n");
            printf("2 - Best Fit\n");
            printf("3 - Worst Fit\n");
            scanf("%u", &code);
            
            if (code2 == 2) {
                printf("Metodo escolhido: First Fit\n\n");
                level_list *resp = first_fit(position_level, level_height, &num_recipients, num_levels, H);
            }
            else if (code2 == 1) {
                printf("Metodo escolhido: Best Fit\n\n");
                level_list *resp = best_fit(position_level, level_height, &num_recipients, num_levels, H);
            }
            else if (code2 == 3) {
                printf("Metodo escolhido: Worst Fit\n\n");
                level_list *resp = worst_fit(position_level, level_height, &num_recipients, num_levels, H);
            }
            else {
                printf("Opcao invalida, digite novamente\n");
                goto line804;
            }

            t = clock() - t; 

            print_tab(code1, code2);
            printf("Quantidade de recipientes: %d\n", num_recipients);
            printf("Tempo de execucao: %f\n", ((double)t)/((CLOCKS_PER_SEC)));

            //armazena a quantidade de recipientes com i niveis na posicao i;
            int recip_with_levels[55];

            for (i = 0; i < 55; i++) {
                recip_with_levels[i]=0;
            }

            for (i = 0; i < num_recipients; i++) {
                recip_with_levels[levels_in_recip[i]]++;
                levels_in_recip[i]=0;
            }

            for (i = 1; i < 51; i++) {
                printf("Recipientes com %d itens: %d\n", i, recip_with_levels[i]);
            }
        }
    }
    return 0;
}



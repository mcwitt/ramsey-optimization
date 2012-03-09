#define TABU_MAXSIZE 1000
#define TABU_HASHLEN 1000

typedef struct
{
    int size;
    int index;
    tabu_elem_t (*hasht)[TABU_HASHLEN];
    tabu_elem_t hist[TABU_MAXSIZE];
} tabu_t;

void tabu_init(tabu_t *tl, int size);   /* call to initialize */

/* add element to history list */
void tabu_add(tabu_t *tl, tabu_elem_t *el);

/* check whether element is in the history list */
int tabu_look(tabu_t *tl, tabu_elem_t *el);

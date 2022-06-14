/*
       ____                _                     ______                                          _____                          _
      / / /     /\        | |                _  |  ____|                                        / ____|                        (_)
     / / /     /  \  _   _| |_ ___  _ __ ___(_) | |__ _ __ __ _ _ __   ___ ___  ___  ___ ___   | |  __ _ __ __ _ _ __   ___ _____  ___
    / / /     / /\ \| | | | __/ _ \| '__/ _ \   |  __| '__/ _` | '_ \ / __/ _ \/ __|/ __/ _ \  | | |_ | '__/ _` | '_ \ / _ \_  / |/ _ \
   / / /     / ____ \ |_| | || (_) | | |  __/_  | |  | | | (_| | | | | (_|  __/\__ \ (_| (_) | | |__| | | | (_| | | | | (_) / /| | (_) |
  /_/_/     /_/    \_\__,_|\__\___/|_|  \___(_) |_|  |_|  \__,_|_| |_|\___\___||___/\___\___/   \_____|_|  \__,_|_| |_|\___/___|_|\___/
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
 * SISTEMA: 
 *      Struttura dati
 * CAMPI:
 *      @a [intero]: rappresenta il coefficiente dell'incognita
 *      @b [intero]: termine noto
 *      @mod [intero]: modulo 
 *      @alpha [intero]: coefficiente di Bézout
 *      @beta [intero]: coefficiente di Bézout
 * 
 * EQUAZIONE: 
 *      ax ≡ b (mod n)
 * 
 * Con a, b appartenenti a Z e d = MCD(a, b).
 * Una coppia di interi(α, β) si dice una coppia di
 * coefficienti di Bézout di a e b se 
 * d = α*a + β*b.
 * 
*/

typedef struct Sistema{

    int a;
    int b;
    int mod;
    int alpha;
    int beta;

}Sistema;


/*
 * INPUT: 
 *      Funzione generica di input di un intero
 * PARAMETRI: 
 *      @msg [const char*]: rappresenta il messaggio da stampare a video 
 * RETURN: 
 *      @var [intero]: intero preso in input
*/

int input(const char *msg){

    int var;
    printf("%s", msg);
    scanf("%d", &var); 

    return var;
}


/*
 * CREASISTEMA: 
 *      Funzione per allocare un array di 
 *      strutture di tipo sistema, richiama "input" prendendo
 *      in input ogni equazione del sistema.
 * PARAMETRI: 
 *      @dim [intero]: numero di equazioni del sistema 
 * RETURN: 
 *      @*sis [puntatore a struttura]: sistema di equazioni
*/

Sistema *creaSistema(int dim){

    Sistema *sis = (Sistema *) malloc(sizeof(Sistema)*dim);

    for(int i = 0; i < dim; i++){

        printf("%d° equazione\n", i + 1);
        sis[i].a = input("Inserisci il coefficiente dell'incognita: ");
        sis[i].b = input("Inserisci il termine noto: ");
        sis[i].mod = input("Inserisci il modulo: ");
    }

    return sis;
}


/*
 * MCD: 
 *      Funzione ricorsiva per calcolare il MCD(a, b)
 * PARAMETRI: 
 *      @a [intero]
 *      @b [intero]
 * RETURN: 
 *      massimo comun divisore tra a e b
*/

int MCD(int a, int b) {

    if (b != 0)
        return MCD(b, a % b);
    else
        return a;
}


/*
 * MODULICOPRIMI: 
 *      Funzione che controlla la compatibilità del sistema.
 * 
 *      {x ≡ b1 (mod n1)
 *      {x ≡ b2 (mod n2)
 *      {...
 *      {x ≡ bt (mod nt)
 * 
 *      Un sistema di equazioni congruenziali ammette sempre soluzioni se 
 *      i moduli sono 2 a 2 coprimi, inoltre se una C appartenente a Z
 *      è una soluzione del sistema, allora le soluzioni del sistema 
 *      saranno tutti e soli gli interi in
 *      [C] MOD n1 * n2 * ... * nt.
 *      Se i moduli non sono 2 a 2 coprimi il programma termina,
 *      altrimenti per il TEOREMA 2 l'equazione congruenziale
 *      ax ≡ b (mod n) è equivalente ad a/d x ≡ b/d (mod n/d),
 *      dove d = MCD(a, n)
 * PARAMETRI: 
 *      @*sis [puntatore a struttura]: sistema di equazioni
 *      @dim [intero]: numero di equazioni del sistema 

*/

void moduliCoprimi(Sistema *sis, int dim){

    int mcd;

    for(int i = 0; i < dim; i++){

        mcd = MCD(sis[i].a, sis[i].mod);
        sis[i].a /= mcd;
        sis[i].b /= mcd;
        sis[i].mod /= mcd;
    }

    for(int i = 0; i < dim; i++){

        for(int j = i + 1; j < dim; j++){

            if(MCD(sis[i].mod, sis[j].mod) != 1){

                puts("I moduli non sono 2 a 2 coprimi, il sistema è incompatibile ");
                exit(EXIT_FAILURE);
            }
        }
    }
}


/*
 * BEZOUT: 
 *      Funzione che trova i coefficienti di Bézout (α, β) tali che 
 *      d = a*α + β*b, dove d = MCD(a, b).
 *      La funzione sfrutta l'algoritmo della divisione euclidea
 * PARAMETRI: 
 *      @a [intero]
 *      @b [intero]
 *      @alpha [intero]: primo coefficiente di Bézout
 *      @beta [intero]: secondo coefficiente di Bézout
*/

void bezout(int a, int b, int *alpha, int *beta){

    int x = 0, y = 1, tmp, quoziente, resto;

    *alpha = 1;
    *beta = 0;

    while (b != 0){

        quoziente = a / b;
        resto = a % b;

        a = b;
        b = resto;

        tmp = x;
        x = *alpha - quoziente * x;
        *alpha = tmp;

        tmp = y;
        y = *beta - quoziente * y;
        *beta = tmp;            
    }
}


/*
 * CLASSEDIEQUIVALENZA: 
 *      Funzione che restitusce la
 *      classe di equivalenza della soluzione C 
 * PARAMETRI: 
 *      @alpha [intero]: primo coefficiente di Bézout
 *      @b [intero]: termine noto
 *      @mod [intero]: modulo equazione
 * RETURN: 
 *      restituisce la classe di equivalenza
 *      della soluzione [C] mod n
*/

int classeDiEquivalenza(int alpha, int b, int mod){

    int num1 = abs(alpha * b);
    int num2 = abs(mod);
    int resto = 0;

    for(int i = 1; i <= num1; i++){

        if(num2 * i > num1){

            resto = num1 - num2 * (i - 1);
            break;
        }
    }

    if(alpha * b < 0)
        resto = mod - resto;

    return resto;
}


/*
 * RIMUOVIA: 
 *      Funzione che trasforma le equazioni del 
 *      sistema da ax ≡ bt (mod nt) a x ≡ ct (mod nt)
 *      dove C è la soluzione dell'equazione
 * PARAMETRI: 
 *      @*sis [puntatore a struttura]: contiene le equazioni del sistema
 *      @dim [intero]: numero equazioni del sistema
*/

void rimuoviA(Sistema *sis, int dim){

    for(int i = 0; i < dim; i++){

        bezout(sis[i].a, sis[i].mod, &sis[i].alpha, &sis[i].beta);
        sis[i].a = 1;
        sis[i].b = classeDiEquivalenza(sis[i].alpha, sis[i].b, sis[i].mod);

    }
}


/*
 * RISOLVISISTEMA: 
 *      Funzione che risolve il sistema di equazioni
 *      Considerando il sistema: 
 * 
 *          {x ≡ b1 (mod n1)
 *          {x ≡ b2 (mod n2)
 *      
 *      Con MCD(n1, n2) = 1.
 * 
 *      Le soluzioni della prima equazione saranno 
 *      b1 + n1 * k, k appartenente a Z.
 *      Adesso bisogna trovare tra tutti i numeri del tipo
 *      b1 + n1 * k almeno uno che soddisfi anche la seconda
 *      equazione, che quindi sarà del tipo
 *      
 *          b1 + n1 * k ≡ b2 (mod n2)
 *      
 *      Da essa si ricava una nuova equazione congruenziale nell'
 *      incognita k.
 *      Adesso, per la proprietà riflessiva, -b1 ≡ -b1 (mod n2), 
 *      sommando queste 2 congruenze modulo n2 si ha se
 *      n1 * k ≡ b2 - b1 (mod n2).
 *      Adesso per l'ipotesi fatta in precedenza MCD(n1, n2) = 1
 *      1 | b2 - b1, e quindi per il TEOREMA 1
 *      l'equazione congruenziale ax ≡ b (mod n) 
 *      ammette soluzioni <==> MCD(a, n) | b.
 *      Quindi una qualunque soluzione di questa equazione
 *      è un valore di k tale che b1 + n1 * k è soluzione del sistema.
 *      
 * PARAMETRI: 
 *      @*sis [puntatore a struttura]: contiene le equazioni del sistema
 *      @dim [intero]: numero equazioni del sistema
*/

void risolviSistema(Sistema *sis, int dim){

    if(dim > 1){

        int k, classe, modulo = sis[0].mod;
        sis[1].b = sis[1].b - sis[0].b;

        for(int i = 1; i < dim; i++){

            sis[i].a = modulo;

            if(i > 1)
                sis[i].b = sis[i].b - classe;

            bezout(sis[i].a, sis[i].mod, &sis[i].alpha, &sis[i].beta);
            k = classeDiEquivalenza(sis[i].alpha, sis[i].b, sis[i].mod);

            if(i > 1)
                classe = sis[i].a * k + classe;
            else
                classe = sis[i - 1].b + sis[i].a * k;

            modulo *= sis[i].mod;

            printf("\nLe soluzioni della %d equazione saranno [%d] MOD %d e del tipo %d + %dk\n", i, classe, modulo, classe, modulo);
        }
    }
    else{

        bezout(sis[0].a, sis[0].mod, &sis[0].alpha, &sis[0].beta);
        int classe =  classeDiEquivalenza(sis[0].alpha, sis[0].b, sis[0].mod);
        printf("Le soluzioni dell'equazione %dx ≡ %d (mod %d) saranno [%d] MOD: %d e del tipo %d + %dk\n", sis[0].a, sis[0].b, sis[0].mod, classe, sis[0].mod, classe, sis[0].mod);
    }

    
}


/*
 * ESEMPIO:
 *  
 *      {36x ≡ 48 (mod 84) => Dividendo per il MCD(36, 84) => 3x ≡ 4 (mod 7)
 *      {14x ≡ 20 (mod 22) => Dividendo per il MCD(14, 22) => 7x ≡ 10 (mod 11)
 *      {21x ≡ 12 (mod 15) => Dividendo per il MCD(21, 15) => 7x ≡ 4 (mod 5)
 * 
 *      Trovando le soluzioni della prima equazione il sistema può essere riscritto come:
 *  
 *      {x ≡ 6 (mod 7)
 *      {x ≡ 3 (mod 11)
 *      {x ≡ 2 (mod 5)
 * 
 *      Il sistema è compatibile, perchè i moduli sono 2 a 2 coprimi 
 *      Le soluzioni saranno nella forma [C] MOD 385
 * 
 *      La prima equazione ha soluzine del tipo 6 + 7 * k, con k in Z,
 *      adesso deve aversi 6 + 7 * k ≡ 3 (mod 11) => 7 * k ≡ -3 (mod 11).
 *      L'equazione ha soluzione per k = 9, sostituendo k si avrà
 *      6 + 7 * 9 = 69, quindi le soluzioni delle prime due equazioni sono 2
 *      [69] MOD 77, del tipo 69 + 77 * k, con k in Z.
 * 
 *      Per la terza equazione deve aversi 69 + 77 * k ≡ 2 (mod 5)
 *      => 77 * k ≡ -67 (mod 5), l'equazione ha soluzione per k = 4.
 *      Sostituendo k si avrà 69 + 77 * 4 = 377.
 * 
 *      Le soluzioni del sistema, dunque, saranno [377] MOD 385
 *      e del tipo 377 + 385 * k, con k in Z.
 *      
 * 
*/


int main(void){

    while(1){

        puts("---------------- RISOLUTORE SISTEMI DI EQUAZIONI CONGRUENZIALI ----------------");
        int dim;
        printf("Dimensione sistema: ");
        scanf("%d", &dim);

        Sistema *sis = creaSistema(dim);

        moduliCoprimi(sis, dim);
        rimuoviA(sis, dim);
        risolviSistema(sis, dim);

        free(sis);
    }
    return 0;
}
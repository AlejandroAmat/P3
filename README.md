PAV - P3: estimación de pitch
=============================

Esta práctica se distribuye a través del repositorio GitHub [Práctica 3](https://github.com/albino-pav/P3).
Siga las instrucciones de la [Práctica 2](https://github.com/albino-pav/P2) para realizar un `fork` de la
misma y distribuir copias locales (*clones*) del mismo a los distintos integrantes del grupo de prácticas.

Recuerde realizar el *pull request* al repositorio original una vez completada la práctica.

Ejercicios básicos
------------------

- Complete el código de los ficheros necesarios para realizar la estimación de pitch usando el programa
  `get_pitch`.

   * Complete el cálculo de la autocorrelación e inserte a continuación el código correspondiente.

      ```c
      void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {
        for (unsigned int l = 0; l < r.size(); ++l) {
          /// \TODO Compute the autocorrelation r[l]
          r[l] =0;
          for(unsigned int n = l; n < x.size(); n++){
              r[l] += x[n]*x[n-l];
          }
          r[l] /= x.size(); //Divide by N
        }

        if (r[0] == 0.0F) //to avoid log() and divide zero 
          r[0] = 1e-10; 
      }
      ```




   * Inserte una gŕafica donde, en un *subplot*, se vea con claridad la señal temporal de un segmento de
     unos 30 ms de un fonema sonoro y su periodo de pitch; y, en otro *subplot*, se vea con claridad la
	 autocorrelación de la señal y la posición del primer máximo secundario.

        Para obtener segmento vocálico junto con la aurtocorrelación que hemos calculado escribimos en 2 ficheros la señal y su autocorrelación.

        ```c
        if(segments == 4){ //4 para que no sea en transición entre sordo y sonoro
          FILE *f_x = fopen("res_x.txt", "w+");
          FILE *f_r = fopen("rex_r.txt", "w+");
          for(unsigned int i =0; i<x.size(); i++){
            fprintf(f_x , "%f \n", x[i]);
            
          }
          for(unsigned int i =0; i<r.size(); i++){
            fprintf(f_r, "%f \n", r[i]);
            
          }
          
          fclose(f_r);
          fclose(f_x);
        }
        ```

        Ya teniendo estos dos ficheros y mediante el código de "autocorr.py" :
        
        ```c
        import numpy as np
        import matplotlib.pyplot as plt

        x = np.loadtxt("res_x.txt", dtype=float)
        r = np.loadtxt("rex_r.txt", dtype=float)

        fm = 20000
        time = np.arange(0,len(x)).astype(float)
        time = time/fm

        muestras = np.arange(0,len(r)).astype(int)



        plt.subplot(211)
        plt.title("Segmento de Señal Sonoro")
        plt.plot(time, x, linewidth =0.5)
        plt.xlabel('Tiempo (s)')
        plt.ylabel('Amplitud')
        plt.grid(True)
        plt.subplot(212)
        plt.xlabel('Tiempo (s)')
        plt.ylabel('Autocorrelación Calculada')
        plt.plot(muestras, r, linewidth =0.5)
        plt.grid(True)
        plt.show()
        }
        ```


      ![Alt text](./img/Autocorr_X.jpg?raw=true "Optional Title")

	NOTA: es más que probable que tenga que usar Python, Octave/MATLAB u otro programa semejante para
	 hacerlo. Se valorará la utilización de la biblioteca matplotlib de Python.

   * Determine el mejor candidato para el periodo de pitch localizando el primer máximo secundario de la
     autocorrelación. Inserte a continuación el código correspondiente.
     
     Vemos por la gráfica que aproximadamente el segundo máximo está en 75. Con lo que 20000/75 = 266 Hz.
     El código:

      ```c
      for (iR=iRMax = r.begin()+npitch_min; iR< r.begin() + npitch_max ;iR++ ){
        if(*iR > *iRMax) iRMax = iR;
      }

      unsigned int lag = iRMax - r.begin();
      //////********************************\\\\\\\
      
      return (float) samplingFreq/(float) lag; ///Pitch
      ```
   * Implemente la regla de decisión sonoro o sordo e inserte el código correspondiente.

      Para la regla de decisión hemos tenido en cuenta los 3 parámetros ya otorgados (nr1, nrmax, pot) y además hemos calculado el ZCR para ser aún más precisos.
      El ZCR:
      ```c
      double zcr =0;
          for(unsigned int i =1; i < x.size();i++){
          if((x[i]>0 && x[i-1]<0) || (x[i]<0 && x[i-1]>0) ) //falta considerar caso 0
            zcr++;
          }
      ```
      Tras esto , nuestro decisor consiste de 4 Umbrales correspondientes a cada una de estas variables. En un principio hemos Hardcodeado los umbrales pero más adelante se explica la optimizacióin.
      ```c
        bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, float zcr) const {
        /// \TODO Implement a rule to decide whether the sound is voiced or not.
        /// * You can use the standard features (pot, r1norm, rmaxnorm),
        ///   or compute and use other ones.
        numero++;
        cout<<"("<<numero+1<<") pot:"<<pot<<", zcr:"<<zcr<<", r1norm:"<<r1norm<<", maxnorm:"<<rmaxnorm<< endl;
        bool unvoiced = true;
        if(r1norm>0.72668 && zcr <149 && pot>-45 &&  rmaxnorm >0.34) unvoiced=false;;

        return unvoiced;
        }
      ```
   * Puede serle útil seguir las instrucciones contenidas en el documento adjunto `código.pdf`.

- Una vez completados los puntos anteriores, dispondrá de una primera versión del estimador de pitch. El 
  resto del trabajo consiste, básicamente, en obtener las mejores prestaciones posibles con él.

  * Utilice el programa `wavesurfer` para analizar las condiciones apropiadas para determinar si un
    segmento es sonoro o sordo. 
	
	  - Inserte una gráfica con la estimación de pitch incorporada a `wavesurfer` y, junto a ella, los 
	    principales candidatos para determinar la sonoridad de la voz: el nivel de potencia de la señal
		(r[0]), la autocorrelación normalizada de uno (r1norm = r[1] / r[0]) y el valor de la
		autocorrelación en su máximo secundario (rmaxnorm = r[lag] / r[0]).

		Puede considerar, también, la conveniencia de usar la tasa de cruces por cero.


        ```c
        FILE *r1 = fopen("r1.txt", "w+");
        FILE *rmax = fopen("rmax.txt", "w+");
        FILE *zcrf = fopen("zcr.txt", "w+");  
        FILE *potf = fopen("pot.txt", "w+");
        fprintf(r1 , "%f \n", r[1]/r[0]);
        fprintf(rmax , "%f \n", r[lag]/r[0]);
        fprintf(zcrf , "%f \n", zcr);
        fprintf(potf , "%f \n", pot);
        ```

        ![Alt text](./img/atribs.jpg?raw=true "Optional Title")

        1) R[1]/R[0]

        2) Potencia

        3) ZCR

        4) Rmax/R[0]

        5) Estimación Pitch (Mencionar que está aplicado el filtro de mediana por eso es tan regular)

	    Recuerde configurar los paneles de datos para que el desplazamiento de ventana sea el adecuado, que
		en esta práctica es de 15 ms.

      - Use el estimador de pitch implementado en el programa `wavesurfer` en una señal de prueba y compare
	    su resultado con el obtenido por la mejor versión de su propio sistema.  Inserte una gráfica
		ilustrativa del resultado de ambos estimadores.

    (I) Wavesurfer:

    ![Alt text](./img/comparativa.jpg?raw=true "Optional Title")

    (II) Python(MatPlotLib):

    ![Alt text](./img/com_pyy.png?raw=true "Optional Title")

    El código es muy parecido al de la antigua gráfica con Pyhton:

      ```c
      import numpy as np
      import matplotlib.pyplot as plt

      x = np.loadtxt("prueba.f0", dtype=float)
      x_ideal = np.loadtxt("prueba.f0ref", dtype=float)


      time = np.arange(0,len(x)).astype(float)
      time = time*0.015




      plt.subplot(211)
      plt.title("Comparativa Estimación Ideal")
      plt.plot(time, x_ideal, linewidth =0.5)
      plt.xlabel('Tiempo (s)')
      plt.ylabel('Frec(Hz)')
      plt.grid(True)
      plt.subplot(212)
      plt.xlabel('Tiempo(s)')
      plt.ylabel('Frec(Hz)')
      plt.plot(time, x, linewidth =0.5)
      plt.grid(True)
      plt.show()
      ```

    Como se observa, son muy parecidas. Nuestra estimación resulta mucho menos irregular dado que hemos aplicado un filtro de mediana (de 3 coeficientes) que posteriormente comentaremos. Por lo demás, son muy parecidas. Si ejecutamos el fichero de evalución obtendremos:
    
    
    ![Alt text](./img/f_score_pr.jpg?raw=true "Optional Title")

    Cabe mencionar que esta puntuación se ha obtenido optimizando la base de datos y no este fichero en concreto. Con lo que no es la puntuación más alta que se podía haber obtenido.
     
	Aunque puede usar el propio Wavesurfer para obtener la representación, se valorará
	el uso de alternativas de mayor calidad (particularmente Python).
  
  * Optimice los parámetros de su sistema de estimación de pitch e inserte una tabla con las tasas de error
    y el *score* TOTAL proporcionados por `pitch_evaluate` en la evaluación de la base de datos 
	`pitch_db/train`..

Ejercicios de ampliación
------------------------

- Usando la librería `docopt_cpp`, modifique el fichero `get_pitch.cpp` para incorporar los parámetros del
  estimador a los argumentos de la línea de comandos.


   
  
  Esta técnica le resultará especialmente útil para optimizar los parámetros del estimador. Recuerde que
  una parte importante de la evaluación recaerá en el resultado obtenido en la estimación de pitch en la
  base de datos.

  * Inserte un *pantallazo* en el que se vea el mensaje de ayuda del programa y un ejemplo de utilización
    con los argumentos añadidos.


    ![Alt text](./img/help.png?raw=true "Optional Title")

    Código:

      ```c
      int main(int argc, const char *argv[]) {
      /// \TODO 
      ///  Modify the program syntax and the call to **docopt()** in order to
      ///  add options and arguments to the program.
      std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
            {argv + 1, argv + argc},	// array of arguments, without the program name
            true,    // show help if requested
            "2.0");  // version string

      std::string input_wav = args["<input-wav>"].asString();
      std::string output_txt = args["<output-txt>"].asString();
      float umbralPot=stof(args["--umbralPot"].asString());
      float umbralZCR=stof(args["--umbralZCR"].asString());
      float umbralR1=stof(args["--umbralR1"].asString());
      float umbralRmax=stof(args["--umbralRmax"].asString());
      ```
- Implemente las técnicas que considere oportunas para optimizar las prestaciones del sistema de estimación
  de pitch.

  Entre las posibles mejoras, puede escoger una o más de las siguientes:

  * Técnicas de preprocesado: filtrado paso bajo, diezmado, *center clipping*, etc.
  * Técnicas de postprocesado: filtro de mediana, *dynamic time warping*, etc.
  * Métodos alternativos a la autocorrelación: procesado cepstral, *average magnitude difference function*
    (AMDF), etc.
  * Optimización **demostrable** de los parámetros que gobiernan el estimador, en concreto, de los que
    gobiernan la decisión sonoro/sordo.
  * Cualquier otra técnica que se le pueda ocurrir o encuentre en la literatura.

  Encontrará más información acerca de estas técnicas en las [Transparencias del Curso](https://atenea.upc.edu/pluginfile.php/2908770/mod_resource/content/3/2b_PS%20Techniques.pdf)
  y en [Spoken Language Processing](https://discovery.upc.edu/iii/encore/record/C__Rb1233593?lang=cat).
  También encontrará más información en los anexos del enunciado de esta práctica.

  Incluya, a continuación, una explicación de las técnicas incorporadas al estimador. Se valorará la
  inclusión de gráficas, tablas, código o cualquier otra cosa que ayude a comprender el trabajo realizado.

  También se valorará la realización de un estudio de los parámetros involucrados. Por ejemplo, si se opta
  por implementar el filtro de mediana, se valorará el análisis de los resultados obtenidos en función de
  la longitud del filtro.
   



- Preprocesado
    * Hamming
    
        Se ha utilizado la ventana de Hamming definida por Matlab para una mejor respuesta de la estimación.


      ![Alt text](./img/hamming.png?raw=true "Optional Title")


        El código necesario:


         ```c
          case HAMMING:
          /// \TODO Implement the Hamming window
     
          if(frameLen%2==0){ //if framelen is even
          unsigned int half = frameLen/2;

            for(unsigned int n =0; n<half; n++){ //Hamming Window for first half
            window[n] = 0.5 * (1-cos(2*M_PI*(n+1)/frameLen+1));
            //fprintf(hamming, "%f \n", window[n]);
            }

            unsigned int idx = half-1;
            for(unsigned int n=half; n<frameLen; n++){ //Symmentric window for the second half
              window[n]= window[idx];
              idx--;
              //fprintf(hamming, "%f \n", window[n]);
            }


          }
          else{

              unsigned int half = (frameLen + 1 )/2;

            for(unsigned int n =0; n<half; n++){ //Hamming Window for first half
            window[n] = 0.5 * (1-cos(2*M_PI*(n+1)/frameLen+1));
            //fprintf(hamming, "%f \n", window[n]);
            }

            int idx = half-2;
            for(unsigned int n=half; n<frameLen; n++){ //Symmentric window for the second half
              window[n]= window[idx];
              idx--;
              //fprintf(hamming, "%f \n", window[n]);
            }
      ```


      
      Si comparamos los resultados de la ventana rectangular con la de Hamming con la base de datos obtenemos:

      * Rectangular:

        ![Alt text](./img/rect_fd.png?raw=true "Optional Title")

      * Hamming:

        ![Alt text](./img/hamming_fd.png?raw=true "Optional Title")


    * Center Clipping



        ![Alt text](./img/clipping.jpg?raw=true "Optional Title")


         ```c
        float max_signal = *std::max_element(x.begin(), x.end());
        float Cl = 0.1133*max_signal;
        for(auto it = x.begin(); it<x.end(); it++){
        if(abs(*it)<Cl) *it=1e-10;
        else if(*it > Cl) *it+=1.317*Cl;
        else if(*it < -Cl) *it-=1.3117*Cl;
        }
         ```

        Comparamos el resultado del Clippping con la base de datos: (Considerar que aún no hemos hecho la optimización)

        * Sin_Clipping

          ![Alt text](./img/not_clipping.png?raw=true "Optional Title")

        * Clipping

          ![Alt text](./img/hamming_fd.png?raw=true "Optional Title")


- PostProcesado
    * Filtro de Mediana de 3 coeficientes

      Hemos considerado el Filtro de mediana y probando 2,3 y 5 coeficientes se comprueba que con 3 es mejor.

      * Código 2 coeficientes (Coge el segundo):
        ```c
        f0_filtered.push_back(f0[0]); //We do not take into account f0(0) anf f0(size-1) as we apply a filter that gets f[n-1], f[n] 
        for(unsigned int l =1; l<f0.size()-1; l++){

        for(int k =-1; k<1; k++){
          filter.push_back(f0[l+k])  ;
          }

        sort(filter.begin(), filter.end());
        f0_filtered.push_back(filter[1]);
        filter.clear();

        }
        ```

      * Con 3 coeficientes :
         ```c
         f0_filtered.push_back(f0[0]); //We do not take into account f0(0) anf f0(size-1) as we apply a filter that gets f[n-1], f[n] and f[n+1];
         for(unsigned int l =1; l<f0.size()-1; l++){

         for(int k =-1; k<2; k++){
            filter.push_back(f0[l+k])  ;
            }

         sort(filter.begin(), filter.end());
         f0_filtered.push_back(filter[1]);
         filter.clear();

        }
        f0_filtered.push_back(f0[f0.size()-1]);
        ```
        * Con 5 coeficientes :
         ```c
         f0_filtered.push_back(f0[0]); //We do not take into account f0(0),f0(1), f0(size-2) anf f0(size-1) as we apply a filter that gets f[n-2], f[n-1], f[n] and f[n+1], f[n+2];
         f0_filtered.push_back(f0[1]); 
         for(unsigned int l =1; l<f0.size()-1; l++){

         for(int k =-2; k<3; k++){
            filter.push_back(f0[l+k])  ;
            }

         sort(filter.begin(), filter.end());
         f0_filtered.push_back(filter[2]);
         filter.clear();

        }
        f0_filtered.push_back(f0[f0.size()-2]);
        f0_filtered.push_back(f0[f0.size()-1]);
        ```


        Comprobando los resultados el que mejor nos da es el filtro de mediana con 3 coeficientes. (Ya se ha tenido en cuenta en las capturas de F_score)


        * Comparativa

          * Gráficas
            ![Alt text](./img/comp_median.png?raw=true "Optional Title")

          * F_score
            Con Filtro: 91.23
            
            Sin filtro: 91.06



- Optimización





Evaluación *ciega* del estimador
-------------------------------

Antes de realizar el *pull request* debe asegurarse de que su repositorio contiene los ficheros necesarios
para compilar los programas correctamente ejecutando `make release`.

Con los ejecutables construidos de esta manera, los profesores de la asignatura procederán a evaluar el
estimador con la parte de test de la base de datos (desconocida para los alumnos). Una parte importante de
la nota de la práctica recaerá en el resultado de esta evaluación.

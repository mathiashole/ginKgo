#!/usr/bin/perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;

# Definir variables para las opciones
my $all_flag = 0;
my $per_flag = 0;
my $tiny_flag = 0;
my $huge_flag = 0;

# Especificar opciones permitidas
GetOptions(
    "all" => \$all_flag,
    "per" => \$per_flag,
    "tiny" => \$tiny_flag,
    "huge" => \$huge_flag
);

# Verifica si se proporcionó al menos un archivo FASTA como argumento
if (@ARGV < 1) {
    die("Uso: $0 <archivo1.fasta> [archivo2.fasta ...]\n");
}

# Itera sobre cada archivo FASTA proporcionado como argumento
foreach my $archivo_fasta (@ARGV) {
#    all_process($archivo_fasta);
if($tiny_flag) {
print "$archivo_fasta\tCodon\tAA\tFraction\tFrequency\tNumber\n";
} elsif($huge_flag) {
    my %tabla = tabla_codon_aa();
    my @codones = sort { $tabla{$a} cmp $tabla{$b} } keys %tabla;

    print "codon\t", join("\t", @codones);
}

if($all_flag) {
    my $secuencia = leer_archivo_fasta($archivo_fasta);
    all_process($secuencia, $archivo_fasta);
} elsif ($per_flag) {
    leer_archivo_fasta($archivo_fasta);
}
} #ESTO PUEDE SER SUB 1 FORMATO 1 Y LA OTRA FORMATO 2

sub all_process {
#my $archivo_fasta = $ARGV[0];
my ($secuencia, $archivo_fasta) = @_;

# Llamar tabla de codones
my %tabla_codon_aa = tabla_codon_aa();

# Inicializa contadores de codones y aminoácidos
$secuencia =~ s/\s//g; # Remueve espacios en blanco
$secuencia = uc($secuencia);
my ($contador_codones, $contador_aa, $total_codones) = contar_codones_y_aa($secuencia, \%tabla_codon_aa);

# Almacena los datos en una lista de hash
my @datos = generar_datos($contador_codones, $contador_aa, $total_codones, \%tabla_codon_aa);

# Ordena la lista por aminoácido
@datos = sort { $a->{aa} cmp $b->{aa} } @datos;

# Imprime la tabla ordenada
imprimir_tabla(\@datos, $archivo_fasta);

# Guarda los resultados en un archivo de texto
guardar_resultados(\@datos, $archivo_fasta);
} ### DIVIDIR ENTRE ALL y MAIN

sub per_process {
    my ($encabezado, $secuencia) = @_;

    my %tabla_codon_aa = tabla_codon_aa();

    # Procesa la secuencia
    $secuencia =~ s/\s//g; # Remueve espacios en blanco
    $secuencia = uc($secuencia);
    my ($contador_codones, $contador_aa, $total_codones) = contar_codones_y_aa($secuencia, \%tabla_codon_aa);

    # Almacena los datos en una lista de hash
    my @datos = generar_datos($contador_codones, $contador_aa, $total_codones, \%tabla_codon_aa);

    # Ordena la lista por aminoácido
    @datos = sort { $a->{aa} cmp $b->{aa} } @datos;

    # Imprime la tabla ordenada
    if($tiny_flag) {
        print_tiny_df_per($encabezado, \@datos);
    } elsif ($huge_flag){
        print_huge_pr_df($encabezado, \@datos);
    } else {
        print_tiny_df_per($encabezado, \@datos);
    }
}

# Subrutina para leer el archivo FASTA
sub leer_archivo_fasta {
    my ($archivo_fasta) = @_;
    open my $fh, '<', $archivo_fasta or die("No se pudo abrir el archivo FASTA: $!\n");
    my $secuencia = "";
    my $encabezado_actual = "";
    while (my $linea = <$fh>) {
        chomp($linea);
        if ($linea =~ /^>/) {
            if ($all_flag) {
                next;
            } elsif ($per_flag) {
                if ($secuencia) {
                    per_process($encabezado_actual, $secuencia);
                    $secuencia = "";
                }
                $encabezado_actual = $linea;
            } else {
                next;
            }
        } else {
            $secuencia .= $linea;
        }
    }
    close $fh;
    return $secuencia;
}

# Subrutina para contar los codones y aminoácidos
sub contar_codones_y_aa {
    my ($secuencia, $tabla_codon_aa) = @_;
    my %contador_codones;
    my %contador_aa;
    my $total_codones = length($secuencia) / 3;
    
    for (my $i = 0; $i < length($secuencia) - 2; $i += 3) {
        my $codon = substr($secuencia, $i, 3);
        if (length($codon) == 3) {
            my $aa = $tabla_codon_aa->{$codon}; # 'X' para codones desconocidos
            $contador_codones{$codon}++;
            $contador_aa{$aa}++;
        }
    }
    
    return (\%contador_codones, \%contador_aa, $total_codones);
}

# Subrutina para generar los datos en formato de lista de hash
sub generar_datos {
    my ($contador_codones, $contador_aa, $total_codones, $tabla_codon_aa, $puntuacion_codones) = @_;
    my @datos;
    
    for my $codon (sort keys %$tabla_codon_aa) {
        my $aa = $tabla_codon_aa->{$codon};
        my $fraction = ($contador_codones->{$codon} // 0) / ($contador_aa->{$aa} // 0.001); # Usar 0.001 para evitar división por cero
        my $frequency = ($contador_codones->{$codon} // 0) / ($total_codones) * 100 // 0;
        my $number = $contador_codones->{$codon} // 0;
        # my $observada = $contador_codones->{$codon} // 0; # 
        # my $esperada = ($contador_aa->{$aa} // 0.001) / (scalar(keys %$tabla_codon_aa)); # 
        # my $puntuacion = $observada / $esperada; #
        push @datos, {
            codon => $codon,
            aa => $aa,
            fraction => $fraction,
            frequency => $frequency,
            number => $number,
            # puntuacion => $puntuacion,
            # observada => $observada,
            # esperada => $esperada,
        };
    }
    
    return @datos;
}

# Subrutina para imprimir la tabla
sub imprimir_tabla {
    my ($datos, $archivo_fasta) = @_;
    #print "Codon\tAA\tFraction\tFrequency\tNumber\n";
    for my $dato (@$datos) {
        printf("%s\t%s\t%s\t%.3f\t%.3f\t%d\n", $archivo_fasta, $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number});
        #printf("%s\t%s\t%.3f\t%.3f\t%d\t%.3f\t%.3f\/%.3f\n", $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number}, $dato->{puntuacion}, $dato->{observada},$dato->{esperada});
    }
}

sub print_tiny_df_per {
    my ($encabezado, $datos) = @_;
    #print "Codon\tAA\tFraction\tFrequency\tNumber\n";
    for my $dato (@$datos) {
        printf("%s\t%s\t%s\t%.3f\t%.3f\t%d\n", $encabezado, $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number});
        #printf("%s\t%s\t%.3f\t%.3f\t%d\t%.3f\t%.3f\/%.3f\n", $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number}, $dato->{puntuacion}, $dato->{observada},$dato->{esperada});
    }
}

# Subrutina para guardar los resultados en un archivo de texto
sub guardar_resultados {
    my ($datos, $archivo_fasta) = @_;
    my $output_nombre_del_fichero = "output_" . $archivo_fasta . ".txt";
    open my $fh, '>', $output_nombre_del_fichero or die("No se pudo abrir el archivo de salida: $!\n");
    print $fh "Fasta\tCodon\tAA\tFraction\tFrequency\tNumber\n";
    for my $dato (@$datos) {
        printf $fh "%s\t%s\t%s\t%.3f\t%.3f\t%d\n", $archivo_fasta, $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number};
        #printf $fh "%s\t%s\t%.3f\t%.3f\t%d\t%.3f\n", $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number}, $dato->{puntuacion};
    }
    close $fh;
}

sub print_huge_pr_df {
    my ($encabezado, $datos) = @_;
    # my %tabla = tabla_codon_aa();
    # my @codones = sort { $tabla{$a} cmp $tabla{$b} } keys %tabla;

    # print "codon\t", join("\t", @codones);
    print "\n$encabezado\t";
    for my $dato (@$datos) {
        printf"%.3f\t", $dato->{frequency};
    }
}

sub tabla_codon_aa {
    my %tabla = (
        'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
        'TGC' => 'C', 'TGT' => 'C',
        'GAC' => 'D', 'GAT' => 'D',
        'GAA' => 'E', 'GAG' => 'E',
        'TTC' => 'F', 'TTT' => 'F',
        'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
        'CAC' => 'H', 'CAT' => 'H',
        'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',
        'AAA' => 'K', 'AAG' => 'K',
        'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',
        'ATG' => 'M',
        'AAC' => 'N', 'AAT' => 'N',
        'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
        'CAA' => 'Q', 'CAG' => 'Q',
        'AGA' => 'R', 'AGG' => 'R', 'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',
        'AGC' => 'S', 'AGT' => 'S', 'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',
        'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
        'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
        'TGG' => 'W',
        'TAC' => 'Y', 'TAT' => 'Y',
        'TAA' => '*', 'TAG' => '*', 'TGA' => '*',
    );
    return %tabla;
}
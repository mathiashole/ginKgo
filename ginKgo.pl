#!/usr/bin/perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;

# Define variables for options
my $all_flag = 0;
my $per_flag = 0;
my $tiny_flag = 0;
my $huge_flag = 0;

# Specify allowed options
GetOptions(
    "all" => \$all_flag,
    "per" => \$per_flag,
    "tiny" => \$tiny_flag,
    "huge" => \$huge_flag
);

# Checks if at least one FASTA file was provided as an argument
if (@ARGV < 1) {
    die("Uso: $0 <archivo1.fasta> [archivo2.fasta ...]\n");
}

# Iterates over each FASTA file given as an argument
foreach my $fasta_file (@ARGV) {
#    all_process($fasta_file);
if($tiny_flag) {
print "$fasta_file\tCodon\tAA\tFraction\tFrequency\tNumber\n";
} elsif($huge_flag) {
    my %tabla = aa_codon_table();
    my @codones = sort { $tabla{$a} cmp $tabla{$b} } keys %tabla;

    print "codon\t", join("\t", @codones);
}

if($all_flag) {
    my $secuencia = read_fasta_file($fasta_file);
    all_process($secuencia, $fasta_file);
} elsif ($per_flag) {
    read_fasta_file($fasta_file);
}
} #ESTO PUEDE SER SUB 1 FORMATO 1 Y LA OTRA FORMATO 2

sub all_process {
#my $fasta_file = $ARGV[0];
my ($secuencia, $fasta_file) = @_;

# Call up codon table
my %aa_codon_table = aa_codon_table();

# Initializes codon and amino acid counters
$secuencia =~ s/\s//g; # Remove white space
$secuencia = uc($secuencia);
my ($contador_codones, $contador_aa, $total_codones) = contar_codones_y_aa($secuencia, \%aa_codon_table);

# Store data in a hash list
my @datos = generar_datos($contador_codones, $contador_aa, $total_codones, \%aa_codon_table);

# Sort the list by amino acid
@datos = sort { $a->{aa} cmp $b->{aa} } @datos;

# Imprime la tabla ordenada
imprimir_tabla(\@datos, $fasta_file);

# Save the results to a text file
guardar_resultados(\@datos, $fasta_file);
} ### DIVIDIR ENTRE ALL y MAIN

sub per_process {
    my ($encabezado, $secuencia) = @_;

    my %aa_codon_table = aa_codon_table();

    # Process the sequence
    $secuencia =~ s/\s//g; # Remove white space
    $secuencia = uc($secuencia);
    my ($contador_codones, $contador_aa, $total_codones) = contar_codones_y_aa($secuencia, \%aa_codon_table);

    # Store data in a hash list
    my @datos = generar_datos($contador_codones, $contador_aa, $total_codones, \%aa_codon_table);

    # Sort the list by amino acid
    @datos = sort { $a->{aa} cmp $b->{aa} } @datos;

    # Print the sorted table
    if($tiny_flag) {
        print_tiny_df_per($encabezado, \@datos);
    } elsif ($huge_flag){
        print_huge_pr_df($encabezado, \@datos);
    } else {
        print_tiny_df_per($encabezado, \@datos);
    }
}

# Subroutine to read FASTA file
sub read_fasta_file {
    my ($fasta_file) = @_;
    open my $fh, '<', $fasta_file or die("No se pudo abrir el archivo FASTA: $!\n");
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

# Subroutine to count codons and amino acids
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

# Subroutine to generate the data in hash list format
sub generar_datos {
    my ($contador_codones, $contador_aa, $total_codones, $tabla_codon_aa, $puntuacion_codones) = @_;
    my @datos;
    
    for my $codon (sort keys %$tabla_codon_aa) {
        my $aa = $tabla_codon_aa->{$codon};
        my $fraction = ($contador_codones->{$codon} // 0) / ($contador_aa->{$aa} // 0.001); # Use 0.001 to avoid division by zero
        my $frequency = ($contador_codones->{$codon} // 0) / ($total_codones) * 100 // 0;
        my $number = $contador_codones->{$codon} // 0;
        # my $observada = $contador_codones->{$codon} // 0; # 
        # my $esperada = ($contador_aa->{$aa} // 0.001) / (scalar(keys %$tabla_codon_aa)); # 
        # my $puntuacion = $observada / $esperada; #
        #my $esperada = $contador_aa->{$aa} // 0.001;

        push @datos, {
            codon => $codon,
            aa => $aa,
            fraction => $fraction,
            frequency => $frequency,
            number => $number,
            # puntuacion => $puntuacion,
            # observada => $observada,
            #esperada => $esperada,
        };
    }
    
    return @datos;
}

# Subroutine to print the table
sub imprimir_tabla {
    my ($datos, $fasta_file) = @_;
    #print "Codon\tAA\tFraction\tFrequency\tNumber\n";
    for my $dato (@$datos) {
        printf("%s\t%s\t%s\t%.3f\t%.3f\t%d\n", $fasta_file, $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number});
        #printf("%s\t%s\t%.3f\t%.3f\t%d\t%.3f\t%.3f\/%.3f\n", $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number}, $dato->{puntuacion}, $dato->{observada},$dato->{esperada});
    }
}

sub print_tiny_df_per {
    my ($encabezado, $datos) = @_;
    #print "Codon\tAA\tFraction\tFrequency\tNumber\n";
    for my $dato (@$datos) {
        printf("%s\t%s\t%s\t%.3f\t%.3f\t%d\n", $encabezado, $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number});
        #printf("%s\t%s\t%.3f\t%.3f\t%d\t%d\n", $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number}, $dato->{esperada});
    }
}

# Subroutine to save results to a text file
sub guardar_resultados {
    my ($datos, $fasta_file) = @_;
    my $output_nombre_del_fichero = "output_" . $fasta_file . ".txt";
    open my $fh, '>', $output_nombre_del_fichero or die("No se pudo abrir el archivo de salida: $!\n");
    print $fh "Fasta\tCodon\tAA\tFraction\tFrequency\tNumber\n";
    for my $dato (@$datos) {
        printf $fh "%s\t%s\t%s\t%.3f\t%.3f\t%d\n", $fasta_file, $dato->{codon}, $dato->{aa}, $dato->{fraction}, $dato->{frequency}, $dato->{number};
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

sub aa_codon_table {
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
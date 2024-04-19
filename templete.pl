@fields = grep { /\S+/ } split /\s+/;
next unless @fields == 2;
next unless \$fields[1] =~ /(\w+).links\.([\w-]+)\.tsv/;
printf qq{%s,%s,%s\n}, \$1, \$2, \$fields[0];

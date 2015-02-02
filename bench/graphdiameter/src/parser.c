#include "parser.h"
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "null.h"
#include "boolean.h"
#include "check_ptr.h"

void error_message(char* descr) {
  fprintf(stderr, "Error in format of description: < %s > SKIPPING\n", descr);
}

int parse_node_descr_to (char *descr, node_t *node) {
  int id;
  int len = strlen(descr);

  // Allocate up to the lenght of the description string for the out and in
  // descriptions.
  char out_s[len], in_s[len];

  int read;

  // FIXME: This is kind of a hack to make the scanf format string
  // work every time. See if we can get rid of it.
  if(descr[len-1] != ' ') {
    char str[len+1];
    strcpy(str, descr);
    str[len] = ' ';
    str[len+1] = '\0';
    read =
        sscanf(str, "%d |%[0123456789 ]|%[0123456789 ]", &id, out_s, in_s);
  } else {
    read =
        sscanf(descr, "%d |%[0123456789 ]|%[0123456789 ]", &id, out_s, in_s);
  }

//  if (read != 3) {
//    error_message(descr);
//    return 5;
//  }

  int n_out = count_numbers(out_s);
  int n_in = count_numbers(in_s);

  if (n_in < 0 || n_out < 0) {
    error_message(descr);
    return 2;
  }

  node_init(node, id, n_out, n_in);

  int inserted = -1;

  inserted = populate_adjacency(out_s, node->out, n_out);
  if (inserted != n_out) {
    error_message(descr);
    node_free(node);
    return 3;
  }

  inserted = populate_adjacency(in_s, node->in, n_in);
  if (inserted != n_in) {
    error_message(descr);
    node_free(node);
    return 4;
  }

  return 0;
}

node_t * parse_node_descr (char *descr) {
  node_t * node = malloc(sizeof(node_t));
  check_ptr(node);

  parse_node_descr_to(descr, node);

  return node;
}

int count_numbers (char *str) {
  int count = 0;
  int i=0;
  int was_digit = FALSE;
  for (; i<strlen(str); ++i) {
    if (isdigit(str[i])) {
      if (!was_digit) {
        ++count;
      }
      was_digit = TRUE;
    } else if (isspace(str[i])) {
      was_digit = FALSE;
    } else {
      return -1;
    }
  }
  return count;
}

int populate_adjacency(char *adj_str, node_id_t *adj, int n) {
  char *str = adj_str;
  node_id_t elem;
  int i = 0;

  while (str != NULL && i<n) {
    elem = strtol(str, &str, 10);
    adj[i] = elem;
    ++i;
  }

  // there are still some elements
  if (str != NULL && strcmp(str, "") != 0) {
    int j = 0;
    for(; j< strlen(str); ++j) {
      if(!isspace(str[j])) {
        fprintf(stderr,
                "%s:%d: ERROR: trying to populate an andjacency list with more "
                "values than the array can host\n", __FILE__, __LINE__);
        return -1;
      }
    }
  }

  return i;
}

/*
 * Reads an entire file into a dinamically allocated buffer.
 */
char * read_file (char * filename) {
  FILE *f = fopen(filename, "rb");
  fseek(f, 0, SEEK_END);
  long pos = ftell(f);
  fseek(f, 0, SEEK_SET);

  char *bytes = malloc(pos);
  // Already checking memory, no need to call check_ptr
  if (bytes == NULL) {
    fprintf(stderr, "ERROR: Not enough memory to load file %s in memory\n",
            filename);
  }
  fread(bytes, pos, 1, f);
  fclose(f);

  return bytes;
}

int count_lines (char *str) {
  int count = 0;
  int len = strlen(str);
  for (int i=0; i < len; ++i) {
    if (str[i] == '\n') {
      ++count;
    }
  }
  if(str[len-1] != '\n') {
    ++count;
  }
  return count;
}

int parse_graph_file (char *filename, node_t **nodes, int *n) {
  char * file_contents = read_file(filename);

  int num_parsed = parse_graph_string(file_contents, nodes, n);

  free(file_contents);
  return num_parsed;
}

int parse_graph_string (char *str, node_t **nodes, int *n) {
  *n = count_lines(str);
  int i = 0;
  int rc = 0;

  if (*n != 0) {
    *nodes = malloc(*n * sizeof(node_t));
    check_ptr(nodes);
  } else {
    fprintf(stderr, "%s:%d: ERROR: trying to allocate 0 bytes.\n",
            __FILE__, __LINE__);
    return -1;
  }

  char *line = strtok(str, "\n");

  if (line != NULL) {
    rc = parse_node_descr_to(line, &((*nodes)[i]));
    if (rc != 0) {
      free(*nodes);
      fprintf(stderr, "Error in parsing node description: rc=%d\n", rc);
      return -1;
    }
    ++i;
  }

  while(line != NULL) {
    line = strtok(NULL, "\n");
    if(line != NULL) {
      rc = parse_node_descr_to(line, &((*nodes)[i]));
      if (rc != 0) {
        free(*nodes);
        fprintf(stderr, "Error in parsing node description: rc=%d\n", rc);
        return -1;
      }
      ++i;
    }
  }

  return i;
}

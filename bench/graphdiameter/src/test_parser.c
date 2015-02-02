#include <check.h>
#include <stdio.h>
#include "parser.h"
#include "graph.h"
#include "null.h"

START_TEST (count_numbers_test) {
  char *str = "1 2 3 4";
  ck_assert_int_eq(count_numbers(str), 4);

  str = "12 3 4 5 2 54  ";
  ck_assert_int_eq(count_numbers(str), 6);

  str = "12 3 4 5 2 string 54  ";
  ck_assert_int_eq(count_numbers(str), -1);
}
END_TEST

START_TEST (populate_adjacency_test) {
  char *str = "1 2 3 4";
  node_id_t adj[4];

  int rc = populate_adjacency(str, adj, 4);
  ck_assert_int_eq(adj[0], 1);
  ck_assert_int_eq(adj[1], 2);
  ck_assert_int_eq(adj[2], 3);
  ck_assert_int_eq(adj[3], 4);
  ck_assert_int_eq(rc, 4);
}
END_TEST

START_TEST (parse_node_descr_test_1) {
  char *descr = "1 | 2 3 4 | 5 6 7";
  node_t *node = parse_node_descr(descr);

  ck_assert (node != NULL);
  
  if (node != NULL) {
    ck_assert_int_eq(node->id, 1);

    ck_assert_int_eq(node->out[0], 2);
    ck_assert_int_eq(node->out[1], 3);
    ck_assert_int_eq(node->out[2], 4);

    ck_assert_int_eq(node->in[0], 5);
    ck_assert_int_eq(node->in[1], 6);
    ck_assert_int_eq(node->in[2], 7);

    node_delete(node);
  }
}
END_TEST

START_TEST (parse_node_descr_test_2) {
  char *descr = "1 | | ";
  node_t *node = parse_node_descr(descr);

  ck_assert (node != NULL);

  if (node != NULL) {
    ck_assert_int_eq(node->id, 1);

    ck_assert_int_eq(node->num_out, 0);
    ck_assert(node->out == NULL);

    ck_assert_int_eq(node->num_in, 0);
    ck_assert(node->in == NULL);

    node_delete(node);
  }
}
END_TEST

START_TEST (parse_node_descr_test_3) {
  char *descr = "1 | 2 | ";
  node_t *node = parse_node_descr(descr);

  ck_assert (node != NULL);

  if (node != NULL) {
    ck_assert_int_eq(node->id, 1);

    ck_assert_int_eq(node->num_out, 1);
    ck_assert_int_eq(node->out[0], 2);

    ck_assert_int_eq(node->num_in, 0);
    ck_assert(node->in == NULL);

    node_delete(node);
  }
}
END_TEST

START_TEST (parse_node_descr_test_4) {
  char *descr = "1 | | 4 ";
  node_t *node = parse_node_descr(descr);

  ck_assert (node != NULL);

  if (node != NULL) {
    ck_assert_int_eq(node->id, 1);

    ck_assert_int_eq(node->num_out, 0);
    ck_assert(node->out == NULL);

    ck_assert_int_eq(node->num_in, 1);
    ck_assert_int_eq(node->in[0], 4);

    node_delete(node);
  }
}
END_TEST

START_TEST (parse_node_descr_to_test_1) {
  char *descr = "2 | |";
  node_t node;

  int rc = parse_node_descr_to(descr, &node);

  ck_assert_int_eq(node.id, 2);
  ck_assert_int_eq(rc, 0);

  ck_assert_int_eq(node.num_out, 0);
  ck_assert(node.out == NULL);

  ck_assert_int_eq(node.num_in, 0);
  ck_assert(node.in == NULL);

  node_free(&node);

}
END_TEST

START_TEST (count_lines_test) {
  int l = count_lines(
        "a\n"
        "sdijb\n"
        "aosjfblk\n"
        "alks\n"
        "lkfn\n"
        "kajs\n"
        "nfk\n"
        "lf");

  ck_assert_int_eq(l, 8);

  l = count_lines(
        "a\n"
        "sdijb\n");

  ck_assert_int_eq(l, 2);

  l = count_lines("a");

  ck_assert_int_eq(l, 1);
}
END_TEST

START_TEST (parse_graph_string_test_1) {
  char str[] =
      "1 | 2 3 | 5\n"
      "2 | 1 3 | 4 6\n";

  node_t *nodes;
  int n;

  int rc = parse_graph_string(str, &nodes, &n);

  ck_assert_int_eq(n, 2);
  ck_assert_int_eq(rc, 2);

  ck_assert_int_eq(nodes[0].id, 1);
  ck_assert_int_eq(nodes[0].num_out, 2);
  ck_assert_int_eq(nodes[0].num_in, 1);
  ck_assert_int_eq(nodes[0].out[0], 2);
  ck_assert_int_eq(nodes[0].out[1], 3);
  ck_assert_int_eq(nodes[0].in[0], 5);

  ck_assert_int_eq(nodes[1].id, 2);
  ck_assert_int_eq(nodes[1].num_out, 2);
  ck_assert_int_eq(nodes[1].num_in, 2);
  ck_assert_int_eq(nodes[1].out[0], 1);
  ck_assert_int_eq(nodes[1].out[1], 3);
  ck_assert_int_eq(nodes[1].in[0], 4);
  ck_assert_int_eq(nodes[1].in[1], 6);

  for (int i=0; i<n; ++i) {
    node_free(&nodes[i]);
  }
  free(nodes);
}
END_TEST

START_TEST (parse_graph_string_test_2) {
  char str[] =
      "1 | 2 3 | 5\n"
      "2 | 1 3 | 4 6";

  node_t *nodes;
  int n;

  int rc = parse_graph_string(str, &nodes, &n);

  ck_assert_int_eq(n, 2);
  ck_assert_int_eq(rc, 2);

  ck_assert_int_eq(nodes[0].id, 1);
  ck_assert_int_eq(nodes[0].num_out, 2);
  ck_assert_int_eq(nodes[0].num_in, 1);
  ck_assert_int_eq(nodes[0].out[0], 2);
  ck_assert_int_eq(nodes[0].out[1], 3);
  ck_assert_int_eq(nodes[0].in[0], 5);

  ck_assert_int_eq(nodes[1].id, 2);
  ck_assert_int_eq(nodes[1].num_out, 2);
  ck_assert_int_eq(nodes[1].num_in, 2);
  ck_assert_int_eq(nodes[1].out[0], 1);
  ck_assert_int_eq(nodes[1].out[1], 3);
  ck_assert_int_eq(nodes[1].in[0], 4);
  ck_assert_int_eq(nodes[1].in[1], 6);

  for (int i=0; i<n; ++i) {
    node_free(&nodes[i]);
  }
  free(nodes);
}
END_TEST

START_TEST (parse_graph_string_test_3) {
  char str[] =
      "0 | 1 2 | 3\n"
      "1 | 2 3 | 4\n"
      "2 | |\n"
      "3 | 2 | 4\n"
      "4 | 5 | 1\n"
      "5 | |";

  node_t *nodes;
  int n;

  int rc = parse_graph_string(str, &nodes, &n);

  ck_assert_int_eq(n, 6);
  ck_assert_int_eq(rc, 6);

  ck_assert_int_eq(nodes[0].id, 0);
  ck_assert_int_eq(nodes[0].num_out, 2);
  ck_assert_int_eq(nodes[0].num_in, 1);
  ck_assert_int_eq(nodes[0].out[0], 1);
  ck_assert_int_eq(nodes[0].out[1], 2);
  ck_assert_int_eq(nodes[0].in[0], 3);

  ck_assert_int_eq(nodes[2].id, 2);
  ck_assert_int_eq(nodes[2].num_out, 0);
  ck_assert_int_eq(nodes[2].num_in, 0);
  ck_assert(nodes[2].out == NULL);
  ck_assert(nodes[2].in == NULL);

  for (int i=0; i<n; ++i) {
    node_free(&nodes[i]);
  }
  free(nodes);
}
END_TEST

START_TEST (parse_graph_string_test_4) {
  char str[] =
      "2 | |\n"
      "5 | |";

  node_t *nodes;
  int n;

  int rc = parse_graph_string(str, &nodes, &n);

  ck_assert_int_eq(n, 2);
  ck_assert_int_eq(rc, 2);

  ck_assert_int_eq(nodes[0].id, 2);
  ck_assert_int_eq(nodes[0].num_out, 0);
  ck_assert_int_eq(nodes[0].num_in, 0);
  ck_assert(nodes[0].out == NULL);
  ck_assert(nodes[0].in == NULL);

  ck_assert_int_eq(nodes[1].id, 5);
  ck_assert_int_eq(nodes[1].num_out, 0);
  ck_assert_int_eq(nodes[1].num_in, 0);
  ck_assert(nodes[1].out == NULL);
  ck_assert(nodes[1].in == NULL);

  for (int i=0; i<n; ++i) {
    node_free(&nodes[i]);
  }
  free(nodes);
}
END_TEST

Suite * parser_suite () {
  Suite *s = suite_create("Parser");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, count_numbers_test);
  tcase_add_test(tc_core, populate_adjacency_test);
  tcase_add_test(tc_core, parse_node_descr_test_1);
  tcase_add_test(tc_core, parse_node_descr_test_2);
  tcase_add_test(tc_core, parse_node_descr_test_3);
  tcase_add_test(tc_core, parse_node_descr_test_4);
  tcase_add_test(tc_core, parse_node_descr_to_test_1);
  tcase_add_test(tc_core, count_lines_test);
  tcase_add_test(tc_core, parse_graph_string_test_1);
  tcase_add_test(tc_core, parse_graph_string_test_2);
  tcase_add_test(tc_core, parse_graph_string_test_3);
  tcase_add_test(tc_core, parse_graph_string_test_4);
  suite_add_tcase(s, tc_core);

  return s;
}

int main() {
  int number_failed = 0;
  Suite *s = parser_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

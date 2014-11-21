add_library(g_local_mu2 SHARED g_local_mu2.c)
set_target_properties(g_local_mu2 PROPERTIES PREFIX "")

add_library(isolin SHARED isolin.c)
set_target_properties(isolin PROPERTIES PREFIX "")

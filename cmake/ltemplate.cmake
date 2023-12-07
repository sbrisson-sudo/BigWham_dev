download_external_project(ltemplate
  URL "https://github.com/GeoEnergyLab-EPFL/LTemplate.git"
  TAG "v${_ltemplate_version}"
  BACKEND GIT
  THIRD_PARTY_SRC_DIR ${_ltemplate_external_dir}
  NO_UPDATE
  )

set(LTEMPLATE_DIR "${_ltemplate_external_dir}/ltemplate" CACHE PATH "Path to LTemplate" FORCE)
mark_as_advanced(LTEMPLATE_DIR)

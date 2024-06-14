;; (
;;  (python-mode . ((mycondaenv . (lambda ()
;;                  (when (string= system-name "mediator.sdsc.edu")
;;                    (setq-local conda-project-env-path (file-name-concat conda-anaconda-home "envs/sa2"))  
;;                    )))))
;;  (python-ts-mode . ((mycondaenv . (lambda ()
;;                  (when (string= system-name "mediator.sdsc.edu")
;;                    (setq-local conda-project-env-path (file-name-concat conda-anaconda-home "envs/sa2"))  
;;                    )))))
;;  (ess-r-mode . ((mycondaenv . (lambda ()
;;                  (when (string= system-name "mediator.sdsc.edu")
;;                    (setq-local conda-project-env-path (file-name-concat conda-anaconda-home "envs/r"))  
;;                    )))))
;;  )

;; (
;;  (python-mode . ((eval . (when (string= system-name "mediator-sdsc.edu")
;;                            (setq conda-project-env-path (file-name-concat conda-anaconda-home "envs/sa2"))))))
;;  )
;; (
;;  (python-mode . ((conda-project-env-path . (file-name-concat conda-anaconda-home "envs/sa2"))))
;;  (python-ts-mode . ((conda-project-env-path . (file-name-concat conda-anaconda-home "envs/sa2"))))
;;  (ess-r-mode . ((conda-project-env-path . (file-name-concat conda-anaconda-home "envs/r"))))
;; )

(
 (python-mode . ((conda-project-env-path . "~/miniforge3/envs/sa2")))
 (python-ts-mode . ((conda-project-env-path . "~/miniforge3/envs/sa2")))
 (ess-r-mode . ((conda-project-env-path . "~/miniforge3/envs/r")))
 )

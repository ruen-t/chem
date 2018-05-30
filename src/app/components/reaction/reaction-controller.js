var host = "http://35.198.250.233/"
var RUN_REACTION_API = host + "upload"
var REACTION_API = host + "resource/reaction"
var GET_OUTPUT_API = host + 'download/'

class ReactionController {
    constructor($scope, $log, $http, $q, $timeout, $window, $mdDialog) {
        console.log("ReactionController")
        var self = this;
        self.simulateQuery = false;
        self.isDisabled = false;
        self.reagentRequired = false;
        self.selectedReactionList = [];
        self.preview = [1, 2, 3, 4, 5, 6]
        self.loopChoice = createLoopChoice();
        // list of `state` value/display objects
        self.states = loadAll();
        self.querySearch = querySearch;
        self.selectedItemChange = selectedItemChange;
        self.searchTextChange = searchTextChange;
        self.sendRequest = sendRequest;
        self.newState = newState;
        self.addReaction = addReaction;
        self.addNewReaction = addNewReaction;
        self.editReaction = editReaction;
        self.removeReaction = removeReaction;
        self.currentSearchedReaction = null;
        self.downloadFile = downloadFile;
        self.outputReady = false;

        $("#sortable").sortable();
        $("#sortable").disableSelection();

        function editReaction(ev) {
            showDialog(ev, "EditReactionDialogController", 'app/components/reaction/dialog/edit-reaction-dialog.html')
        }
        function addNewReaction(ev) {
            showDialog(ev, "DialogController", 'app/components/reaction/dialog/new-reaction-dialog.html')
        }
        function showDialog(ev, controller, html){
            $mdDialog.show({
                controller: controller + " as vm",
                templateUrl: html,
                parent: angular.element(document.body),
                targetEvent: ev,
                clickOutsideToClose: true,
                fullscreen: $scope.customFullscreen // Only for -xs, -sm breakpoints.
            })
        }
        function addReaction() {
            if (self.currentSearchedReaction != null) {
                self.selectedReactionList.push(self.currentSearchedReaction);
                self.currentSearchedReaction = null;
            }

        }
        function removeReaction(id) {
            console.log("remove: " + id)
            for (var i = 0; i < self.selectedReactionList.length; i++) {
                if (self.selectedReactionList[i].id == id) {
                    self.selectedReactionList.splice(i, 1);
                    break;
                }
            }
        }
        function downloadFile() {
            $window.open(GET_OUTPUT_API + self.outputFile, '_blank');
        }
        function sendRequest() {

            var payload = new FormData();
            payload.append("file", $scope.file);
            var reaction = [];
            self.selectedReactionList.forEach(function (item) {
                reaction.push({ reaction: item.id, reagent: item.input_reagent, loop: item.loop })
            });
            payload.append("reaction", JSON.stringify(reaction));

            $http({
                url: RUN_REACTION_API,
                method: 'POST',
                data: payload,
                //assign content-type as undefined, the browser
                //will assign the correct boundary for us
                headers: { 'Content-Type': undefined },
                //prevents serializing payload.  don't do it.
                transformRequest: angular.identity
            }).then(function (response) {
                console.log(response);
                try {
                    var data = response.data;
                    console.log(data)
                    self.outputFile = data.output;
                    self.outputReady = true;
                    if (typeof data.sample != 'undefined') {
                        self.sampleResult = data.sample;
                    }
                } catch (error) {
                    console.log(error);
                }

            }, function (error) {
                console.log(error)
            })

        }
        function loadAll() {
            return new Promise(function (resolve, reject) {
                $http.get(REACTION_API).then(function (response) {
                    self.reactionList = response.data.reaction;
                    self.reactionList.forEach(function (element) {
                        element.display = element.name + " " + element.smart;
                        element.loop = 1;
                    });
                    resolve()
                });
            })
        }
        function querySearch(query) {
            //console.log(self.reactionList)
            if (typeof self.reactionList == 'undefined') {
                var deferred = $q.defer();
                loadAll().then((success) => {
                    var results = query ? self.reactionList.filter(createFilterFor(query)) : self.reactionList;
                    deferred.resolve(results)
                })
                return deferred.promise;
            } else {
                var results = query ? self.reactionList.filter(createFilterFor(query)) : self.reactionList;
                return results

            }



        }
        function searchTextChange(text) {
            $log.info('Text changed to ' + text);
        }

        function selectedItemChange(item) {
            if (item != null) {
                self.currentSearchedReaction = item;
                addReaction();
                //self.selectedItem = null;
                //self.searchText = "";
                $log.info('Item changed to ' + JSON.stringify(item));
            }


        }
        function newState(state) {
            alert("Sorry! You'll need to create a Constitution for " + state + " first!");
        }
        function createFilterFor(query) {
            console.log("create filter");
            var lowercaseQuery = query.toLowerCase();

            return function filterFn(state) {

                return (state.name.includes(lowercaseQuery));
            };
        }
        function createLoopChoice() {
            var choice = [];
            for (var i = 1; i <= 5; i++) {
                choice.push({ display: i, value: i });
            }
            return choice;
        }


    }

}

angular
    .module('app')
    .component('reaction', {
        templateUrl: 'app/components/reaction/index.html',
        controller: ReactionController,
        controllerAs: "vm"
    })
    .directive('chooseFile', function () {
        return {
            link: function (scope, elem, attrs) {
                var button = elem.find('button');
                var input = angular.element(elem[0].querySelector('input#fileInput'));
                button.bind('click', function () {
                    input[0].click();
                });
                input.bind('change', function (e) {
                    scope.$apply(function () {
                        var files = e.target.files;
                        if (files[0]) {
                            scope.file = files[0];
                            scope.fileName = files[0].name;

                        } else {
                            scope.fileName = null;
                        }
                    });
                });
            }
        };
    });